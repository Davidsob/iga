#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "base/SimulationResults.h"
#include "base/StoredVariable.h"
#include "base/StoredVariableOperator.h"
#include "base/Values.h"
#include "base/VariableManager.h"

#include "iga/IntegrationPoints.h"
#include "iga/IgaIO.h"
#include "iga/ParametricMesh.h"
#include "iga/Quadrature.h"
#include "iga/QuadratureMesh.h"
#include "iga/SolidElementMapper.h"

#include "iga/base/LaSystem.h"

#include "iga/constraints/ConstraintManager.h"
#include "iga/constraints/DirichletConstraint.h"

#include "iga/solvers/LinearSolver.h"
#include "iga/solvers/ContinuationSolver.h"
#include "iga/solvers/NonlinearSolver.h"

#include "iga/weakforms/WeakForms.h"
#include "iga/weakforms/WeakFormTags.h"
#include "iga/weakforms/WeakFormManager.h"

#include "splines/Nurbs.h"
#include "splines/utils/Transformations.h"

#include "thermal/StoredTemperatureVariable.h"
#include "thermal/StoredTemperatureRateVariable.h"

#include "utils/MatrixTypes.h"

#include <Eigen/SparseLU>

static const std::string simulation_name = "ht3d";

/*
/////////////////////////////////////////////////////////////////////////////////////
GEOMETRY
/////////////////////////////////////////////////////////////////////////////////////
*/
void circle(NurbsCurve &curve, double radius=1.0)
{
  using namespace vector_ops;
  int p = 2;
  std::vector<double> knot{0.0, 0.0, 0.0, 1, 1, 2, 2, 3, 3, 4, 4, 4};
  algo::normalizeKnot(knot);
  typename NurbsCurve::matrix cpts{
    {radius , 0.0    , 0.0},
    {radius , radius , 0.0},
    {0.0    , radius , 0.0},
    {-radius, radius , 0.0},
    {-radius, 0.0    , 0.0},
    {-radius, -radius, 0.0},
    {0.0    , -radius, 0.0},
    {radius , -radius, 0.0},
    {radius , 0.0    , 0.0},
  };

  static double const fact{1.0/std::sqrt(2.0)};
  decltype(knot) weights(cpts.size(),1);
  for (size_t i = 1; i < weights.size(); i+= 2)
    weights[i] = fact; 
  curve.p = p;
  curve.knot = knot;
  curve.weights = weights;
  curve.Q = cpts;
}

void sweepCurve(NurbsCurve &curve, double radius=1.0)
{
  using namespace vector_ops;
  int p = 2;
  std::vector<double> knot{0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0};
  algo::normalizeKnot(knot);
  typename NurbsCurve::matrix cpts{
    {0.0 , 0.0    , 4.0},
    {0.0 , 0.0    , 2.0},
    {0.0 , 0.0    , 0.0},
    {0.0 , 0.0    , -radius},
    {radius     , 0.0    , -radius},
  };

  // compute weights
  static double const fact{1.0/std::sqrt(2.0)};
  decltype(knot) weights(cpts.size(),1.0);
  for (size_t i = 3; i < weights.size(); i+= 2) weights[i] = fact; 

  // set properties of curve
  curve.p = p;
  curve.knot = knot;
  curve.weights = weights;
  curve.Q = cpts;
}

void annulus(NurbsSurface &surf, double ri, double ro)
{
  using namespace spline_ops;

  NurbsCurve ci, co;
  circle(ci,ri);
  circle(co,ro);

  surf.p = 2;
  surf.q = 1;
  surf.Q = ci.Q;
  surf.Q.insert(surf.Q.end(), co.Q.begin(), co.Q.end());
  surf.weights = ci.weights;
  surf.weights.insert(surf.weights.end(), co.weights.begin(), co.weights.end());

  surf.uknot = co.knot;
  surf.vknot = std::vector<double>({0,0,1,1});
}

void elbow(double ri, double ro, double re, NurbsSolid &solid)
{
  NurbsCurve curve; sweepCurve(curve,re);
  NurbsSurface surf; annulus(surf,ri,ro);

  solid.p = surf.p;
  solid.q = surf.q;
  solid.r = curve.p;

  solid.uknot = surf.uknot;
  solid.vknot = surf.vknot;
  solid.wknot = curve.knot;

  // now the fun part
  auto addPoints = [&solid](double cw, auto const &section)
  {
    solid.Q.insert(solid.Q.end(), section.Q.begin(), section.Q.end());
    for (auto const &w : section.weights){
      solid.weights.push_back(w*cw);
    }
  };

  {
    auto ders = spline_ops::CurveDerivatives(0,1,curve);
    auto section = transform::project(surf,{0,0,1},ders[0], ders[1]);
    addPoints(curve.weights[0],section);
  }

  {
    auto ders = spline_ops::CurveDerivatives(0.25,1,curve);
    auto section = transform::project(surf,{0,0,1},ders[0], ders[1]);
    addPoints(curve.weights[1],section);
  }

  {
    auto ders = spline_ops::CurveDerivatives(0.5,1,curve);
    auto section = transform::project(surf,{0,0,1},ders[0], ders[1]);
    addPoints(curve.weights[2],section);
  }

  {
    auto ders = spline_ops::CurveDerivatives(0.75,1,curve);
    auto section = transform::project(surf,{0,0,1},ders[0], ders[1]);
    addPoints(curve.weights[3],section);
  }

  {
    auto ders = spline_ops::CurveDerivatives(1.0,1,curve);
    auto section = transform::rotate(surf,{0,0,1}, ders[0], ders[1]);
    addPoints(curve.weights[4],section);
  }
}

/*
/////////////////////////////////////////////////////////////////////////////////////
Handle Simulation Results
/////////////////////////////////////////////////////////////////////////////////////
*/

struct SurfaceWriter
{
  SurfaceWriter(size_t x, size_t y)
  : x(x), y(y){};

  template<typename T>
  void operator()(NurbsSurface const &shape, T const &data, std::string const fname) const
  {
    IO::writeSolutionToFile(shape,data,fname,x,y);
  }
  size_t x,y;
};

struct SolidWriter
{
  SolidWriter(size_t x, size_t y, size_t z)
  : x(x), y(y), z(z) {};

  template<typename T>
  void operator()(NurbsSolid const &shape, T const &data, std::string const fname) const
  {
    IO::writeSolutionToFile(shape,data,fname,x,y,z);
  }
  size_t x,y,z;
};

template<typename Shape, typename IgaWriter>
class Results : public SimulationResults
{
public:

  Results(Shape const &shape, IgaWriter const &writer, std::string simName, size_t freq = 1)
    : _shape(shape)
    , _writer(writer)
    ,_simName(simName)
    , _freq(freq)
    , _step(1)
  {}

  ~Results() = default;

  void write() const override
  {
    if (_step == 1 || (_step%_freq == 0))
    {
      auto time = SimulationClock::instance().time();
      std::cout << "* Simulation time = " << time << " [s]" << std::endl;
      writeTime();
      writeStoredVar<StoredTemperatureVariable>("T");
      // writeStoredVar<StoredExtrapolatedElectricFieldVariable<2>>
      // ("Ey", [](StaticVectorR<2> const &x) { return x[1];});
    }
    _step++;
  }

private:

  void writeTime() const
  {
    static std::string const path      =  "output/";
    static std::string const extension = ".txt";

    auto fname = _simName + "_time";
    std::ofstream ofs;
    ofs.precision(8);
    if (_step == 1) {
      ofs.open(path + fname + extension);
    } else {
      ofs.open(path + fname + extension, std::ofstream::out | std::ofstream::app);
    }
    ofs << SimulationClock::instance().time() << std::endl;
    ofs.close();
  }

  template<typename Var, typename Op>
  void writeStoredVar(std::string const &varname, Op const &&op) const
  {
    static std::string const path      =  "output/";
    static std::string const extension = ".txt";

    std::string const fname = _simName + "_" + varname + "_" + std::to_string(_step);
    auto const &geom = GeometricDofManager::instance();
    auto const &var = *(VariableManager::instance().has<Var>());

    std::vector<double> x;
    for (auto id : geom.idsForShape(&_shape)) x.push_back(op((var)[id]));

    _writer(_shape,x,path+fname+extension);
  }

  template<typename Var>
  void writeStoredVar(std::string const &varname) const
  {
    static std::string const path      =  "output/";
    static std::string const extension = ".txt";

    std::string const fname = _simName + "_" + varname + "_" + std::to_string(_step);
    auto const &geom = GeometricDofManager::instance();
    auto const &var = *(VariableManager::instance().has<Var>());

    std::vector<double> x;
    for (auto id : geom.idsForShape(&_shape)) x.push_back((var)[id]);

    _writer(_shape,x,path+fname+extension);
  }

  Shape const &_shape;
  IgaWriter const &_writer;
  std::string const _simName;
  size_t const _freq;
  mutable size_t _step;
};

/*
/////////////////////////////////////////////////////////////////////////////////////
OPERATORS
/////////////////////////////////////////////////////////////////////////////////////
*/

class Conductivity 
{
public:
  using value_t = double;
  Conductivity() = default; 
  ~Conductivity() = default;

  template<typename Point>
  value_t operator()(Point const &p) const
  {
    return value_t(1.0);
  }
};

class GradOperator 
{
public:
  using value_t = DynamicMatrixR;
  GradOperator() = default; 
  ~GradOperator() = default;

  template<typename Point>
  value_t operator()(Point const &p) const
  {
    return p.mapper.grad(p.para);
  }
};

/*
/////////////////////////////////////////////////////////////////////////////////////
Constraint
/////////////////////////////////////////////////////////////////////////////////////
*/

template<typename ValueType>
class AppliedTemperature
  : public DirichletConstraint<ValueType, TemperatureVariable>
{
public:
  template<typename ...Args>
  AppliedTemperature(Args ...args) : DirichletConstraint<ValueType, TemperatureVariable>(args...) {}
  ~AppliedTemperature() = default;
};

static std::string const python{"~/anaconda2/bin/python2.7 "};
static std::string const geometryFile{"elbow_geometry.txt"};

template<typename Shape> void show(Shape const &shape, std::string const &file)
{
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());
}

/*
/////////////////////////////////////////////////////////////////////////////////////
Simulation
/////////////////////////////////////////////////////////////////////////////////////
*/

void sim()
{
  // get the geometry
  NurbsSolid shape; elbow(1.0, 2.0, 3.0, shape);
  // write to file
  // spline_ops::writeToFile(shape,geometryFile,20,2,20);
  // show geometry
  // show(shape,geometryFile);
  {
    // configure the weak forms
    using bdb_form = BdBWeakForm<GradOperator, Conductivity, DynamicMatrixR>;
    // add the weak forms
    auto &wfm = WeakFormManager::instance();
    wfm.add<bdb_form>(GradOperator(), Conductivity());
  } 
  // add variable storage
  {
    auto mgr = GeometricDofManager::instance();
    mgr.addShape(&shape);
    auto &var = VariableManager::instance();
    auto cpts = mgr.ids.size(); // total of unique cpts
    var.add<StoredTemperatureVariable>(cpts);
    var.add<OldStoredTemperatureVariable>(cpts);
    var.add<StoredTemperatureVariableIncrement>(cpts);
    var.add<StoredTemperatureRateVariable>(cpts);
    var.add<OldStoredTemperatureRateVariable>(cpts);
  }

  // set up the boundary conditinos
  NurbsSurface  hot; spline_ops::getIsoSurface(0,2,shape,hot);
  NurbsSurface cold; spline_ops::getIsoSurface(1,2,shape,cold);
  {
    using ramped_T = AppliedTemperature<RampedValue<double>>;
    using fixed_T = AppliedTemperature<ConstantValue<double>>;
    auto &mgr = ConstraintManager::instance();
    mgr.add<ramped_T>(hot,300.0,1.0,200.0);
    mgr.add<fixed_T>(cold,200.0);
  }

  // solve
  {
    using lasystem_t = FeLaSystem<StoredTemperatureVariable>;
    using inner_t = NonlinearSolver<lasystem_t, StoredTemperatureVariable>;
    using solver_t = ContinuationSolver<inner_t>;
    double const tol{1e-8};
    size_t const inner_steps{15};
    size_t const load_increments{5};
    solver_t solver(tol,inner_steps,&shape);
    solver.setStepSize(1.0/load_increments);

    // add results handling
    using results_t = Results<NurbsSolid,SolidWriter>;
    solver.addResultsHandler<results_t>(shape, SolidWriter(12,3,24), simulation_name);

    solver.run();
  }
}

int main(int argc, char **argv)
{
  sim();
  return 0;
}
