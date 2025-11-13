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
#include "iga/ShapeFunctions.h"
#include "iga/ManifoldElementMapper.h"
#include "iga/ParametricMesh.h"
#include "iga/Quadrature.h"
#include "iga/QuadratureMesh.h"

#include "iga/base/LaSystem.h"

#include "iga/boundary_conditions/BoundaryConditionManager.h"

#include "iga/constraints/ConstraintManager.h"
#include "iga/constraints/DirichletConstraint.h"

#include "iga/solvers/LinearSolver.h"
#include "iga/solvers/ContinuationSolver.h"
#include "iga/solvers/NonlinearSolver.h"

#include "iga/weakforms/WeakForms.h"
#include "iga/weakforms/WeakFormTags.h"
#include "iga/weakforms/WeakFormManager.h"

#include "materials/linear_elastic/MembraneMaterialTangent.h"
#include "materials/linear_elastic/TransverseShearMaterialTangent.h"

#include "splines/Nurbs.h"
#include "splines/SplineModifiers.h"

#include "shellmechanics/boundary_conditions/AppliedPressure.h"
#include "shellmechanics/boundary_conditions/AppliedTraction.h"
#include "shellmechanics/boundary_conditions/FixedFixed.h"
#include "shellmechanics/boundary_conditions/NormalSourceBase.h"
#include "shellmechanics/boundary_conditions/DirichletPointConstraint.h"
#include "shellmechanics/boundary_conditions/DisplacementSymmetryConstraint.h"
#include "shellmechanics/boundary_conditions/RotationSymmetryConstraint.h"
#include "shellmechanics/boundary_conditions/SymmetryConstraint.h"
#include "shellmechanics/boundary_conditions/PointLoad.h"
#include "shellmechanics/variables/StoredSixDofDisplacementVariable.h"
#include "shellmechanics/variables/StoredSixDofVelocityVariable.h"
#include "shellmechanics/operators/DisplacementIdOperator.h"
#include "shellmechanics/operators/LinearShellConstraint.h"
#include "shellmechanics/operators/ShellWeakForms.h"

#include "utils/GeometryUtils.h"
#include "utils/MatrixTypes.h"
#include "utils/MatrixTypeTraits.h"

#include <Eigen/SparseLU>

static std::string const python{"~/anaconda2/bin/python2.7 "};
static std::string const geometryFile{"output/elbow_geometry.txt"};
static std::string const simulation_name = "shell3d-beam";

using primary_var_t = StoredSixDofDisplacementVariable;
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
      writeStoredVar<primary_var_t>("U", [](auto const &x) { return x[2]; });
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
    auto const &sids = geom.idsForShape(&_shape);
    std::vector<double> x;

    // add displacement!
    double scale = 1.0;
    Shape deformed = _shape;
    for (size_t i = 0; i < sids.size(); i++)
    {
      auto const id = sids[i];
      DynamicVectorR const ui = var[id].head(3);
      x.push_back(op(ui));
      deformed.Q[i] += scale*convert::to<std::vector<double>>(ui);
    }

      _writer(deformed,x,path+fname+extension);

      {
        DynamicMatrixR uvw = convert::to<DynamicMatrixR>(var.data());
        auto mx = uvw.colwise().maxCoeff();
        auto mn = uvw.colwise().minCoeff();
        std::cout << "\n#### Solution range:" << std::endl;
        std::printf("      u :{%10.4g, %10.4g}\n", mn[0],mx[0]);
        std::printf("      v :{%10.4g, %10.4g}\n", mn[1],mx[1]);
        std::printf("      w :{%10.4g, %10.4g}\n", mn[2],mx[2]);
        std::printf("      a1:{%10.4g, %10.4g}\n", mn[3],mx[3]);
        std::printf("      a2:{%10.4g, %10.4g}\n", mn[4],mx[4]);
        std::printf("      a3:{%10.4g, %10.4g}\n", mn[5],mx[5]);
      }
    // for (auto id : sids x.push_back(op(var[id]));
    }

  Shape const &_shape;
  IgaWriter const &_writer;
  std::string const _simName;
  size_t const _freq;
  mutable size_t _step;
};

/*
/////////////////////////////////////////////////////////////////////////////////////
GEOMETRY
/////////////////////////////////////////////////////////////////////////////////////
*/
void acylinder(double R, double L, NurbsSurface &surface)
{
  surface.p = 2;
  surface.q = 1;
  surface.uknot = {0,0,0,1,1,1};
  surface.vknot = {0,0,1,1};

  double const fact{1.0/std::sqrt(2.0)};

  surface.Q = {{0,0,R},{R,0,R},{R,0,0},
               {0,L,R},{R,L,R},{R,L,0}};
  surface.weights = {1,fact,1,1,fact,1};
}

template<typename Shape> void show(Shape const &shape, std::string const &file)
{
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());
}

template<typename Shape> void showCurve(Shape const &shape, std::string const &file)
{
  std::system(std::string(python + "python/plot_curve.py " + file).c_str());
}

/*
/////////////////////////////////////////////////////////////////////////////////////
Simulation
/////////////////////////////////////////////////////////////////////////////////////
*/

struct NormalSource : public NormalSourceBase
{
  explicit NormalSource(NurbsSurface const &shape)
    : _shape(shape)
  {}

  ~NormalSource() = default;

  StaticVectorR<3> const operator()(double const u, double const v) const override
  {
    auto const dxd1 = spline_ops::SurfaceDerivative(u,v,1,0,_shape);
    auto const dxd2 = spline_ops::SurfaceDerivative(u,v,1,1,_shape);
    auto const n = normalize(cross(dxd1,dxd2));
    return StaticVectorR<3>({n[0],n[1],n[2]});
  }

  NurbsSurface const &_shape;
};

struct EdgeNormalSource : public NormalSourceBase
{
  explicit EdgeNormalSource(StaticVectorR<3> const &n) : n(n) {}
  ~EdgeNormalSource() = default;

  StaticVectorR<3> const operator()(double const u, double const v) const override
  {
    return n;
  }

  StaticVectorR<3> const n;
};

template<typename Plane>
std::vector<BSplinePoint> pointsOnPlane(Plane const &plane)
{
  auto &geom = GeometricDofManager::instance();
  std::vector<BSplinePoint> found;
  for (auto const id : geom.ids)
  {
    if (GeometryUtils::on_plane(convert::to<DynamicVectorR>(geom.ctrlPoints[id]),plane))
      found.push_back(BSplinePoint(geom.ctrlPoints[id]));
  }
  return found;
}

void sim()
{
  // get the geometry
  double const R = 1; double const L = 0.2; double const thickness = 0.001;
  NurbsSurface shape; acylinder(R,L,shape);
  spline_ops::elevate(2,0,shape);
  spline_ops::elevate(1,1,shape);
  spline_ops::midpointRefinement(5,0,shape);
  // spline_ops::midpointRefinement(3,1,shape);
  {
    auto &geom = GeometricDofManager::instance();
    geom.addShape(&shape);
  }
  // // show geometry
  // spline_ops::writeToFile(shape,geometryFile,10,10);
  // show(shape,geometryFile);

  // the material
  double const E = 80e9; double const nu = 0.0;
  MembraneMaterialTangent matl_m(E,nu);
  TransverseShearMaterialTangent matl_t(E,nu);
  {
    // add the weak forms
    using material_form = LinearMaterialStiffness<MembraneMaterialTangent,TransverseShearMaterialTangent>;
    auto &wfm = WeakFormManager::instance();
    wfm.add<material_form>(material_form(matl_m, matl_t, thickness));
    // residual
    using residual_form = MaterialResidual<MembraneMaterialTangent,TransverseShearMaterialTangent>;
    wfm.add<residual_form>(residual_form(matl_m, matl_t, thickness));
    // add the constraint for the normal displacements
    auto &cm = ConstraintManager::instance();
    cm.add<LinearShellConstraint>(shape);
  } 
  // add variable storage
  {
    auto mgr = GeometricDofManager::instance();
    mgr.addShape(&shape);
    auto &var = VariableManager::instance();
    auto cpts = mgr.ids.size(); // total of unique cpts
    var.add<primary_var_t>(cpts);
    var.add<OldStoredSixDofDisplacementVariable>(cpts);
    var.add<StoredSixDofDisplacementVariableIncrement>(cpts);
    var.add<StoredSixDofVelocityVariable>(cpts);
    var.add<OldStoredSixDofVelocityVariable>(cpts);
  }

  // set up symmetry conditions 
  NurbsCurve xsymm,zsymm;
  spline_ops::getIsoCurve(0.0,0,shape,xsymm);
  spline_ops::getIsoCurve(1.0,0,shape,zsymm);
  // spline_ops::writeToFile(ysymm,geometryFile,20);
  // showCurve(ysymm, geometryFile);
  {
    // auto &mgr = ConstraintManager::instance();
    // mgr.add<FixedFixed>(zsymm);
  }

  GeometryUtils::Plane<3> plane({0,0,0},{0,0,-1});
  auto const onz = pointsOnPlane(plane);
  {
    auto &mgr = ConstraintManager::instance();
    for (auto const &p : onz) mgr.add<PrescribedShellDisplacementU>(p);
    for (auto const &p : onz) mgr.add<PrescribedShellDisplacementV>(p);
    for (auto const &p : onz) mgr.add<PrescribedShellDisplacementW>(p);
    for (auto const &p : onz) mgr.add<PrescribedShellDisplacementA1>(p);
    for (auto const &p : onz) mgr.add<PrescribedShellDisplacementA2>(p);
    for (auto const &p : onz) mgr.add<PrescribedShellDisplacementA3>(p);
  }

  {
    using value_t = StaticVectorR<3>;
    auto &mgr = BoundaryConditionManager::instance();
    using load_T = AppliedTraction<ConstantValue<value_t>>;
    // using load_T = AppliedTraction<RampedValue<double>>;
    mgr.add<load_T>(xsymm,value_t({0,0,-5.0}));
  }

  // solve
  {
    using lasystem_t = FeLaSystem<primary_var_t>;
    using solver_t = LinearSolver<lasystem_t,primary_var_t>;
    solver_t solver(&shape);
    solver.run();

    // using inner_t = NonlinearSolver<lasystem_t, primary_var_t>;
    // using solver_t = ContinuationSolver<inner_t>;

    // double const tol{1e-8};
    // size_t const inner_steps{20};
    // size_t const load_increments{5};
    // solver_t solver(tol,inner_steps,&shape);
    // solver.setStepSize(1.0/load_increments);
    // // add results handling
    // using results_t = Results<NurbsSurface,SurfaceWriter>;
    // solver.addResultsHandler<results_t>(shape, SurfaceWriter(10,10), simulation_name);
    // solver.run();
    using results_t = Results<NurbsSurface,SurfaceWriter>;
    results_t results(shape, SurfaceWriter(50,10), simulation_name);
    results.write();
  }



}

int main(int argc, char **argv)
{
  sim();
  return 0;
}
