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

#include "iga/boundary_conditions/BoundaryConditionManager.h"

#include "iga/constraints/ConstraintManager.h"
#include "iga/constraints/DirichletConstraint.h"

#include "iga/solvers/LinearSolver.h"
#include "iga/solvers/ContinuationSolver.h"
#include "iga/solvers/NonlinearSolver.h"

#include "iga/weakforms/WeakForms.h"
#include "iga/weakforms/WeakFormTags.h"
#include "iga/weakforms/WeakFormManager.h"

#include "splines/Nurbs.h"
#include "splines/SplineModifiers.h"
#include "splines/utils/Transformations.h"

#include "solidmechanics/StoredDisplacementVariable.h"
#include "solidmechanics/StoredVelocityVariable.h"

#include "utils/MatrixTypes.h"
#include "utils/MatrixTypeTraits.h"

// #include <Eigen/SparseLU>

static const std::string simulation_name = "nlsm3d";
using primary_var_t = StoredDisplacementVariable<3>;
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

void abeam(double L, double H, double D, NurbsSolid &solid)
{
  solid.p = 1;
  solid.q = 1;
  solid.r = 1;

  solid.uknot = {{0,0,1,1}};
  solid.vknot = {{0,0,1,1}};
  solid.wknot = {{0,0,1,1}};

  solid.Q = {
    {0,0,0}, {L,0,0},
    {0,H,0}, {L,H,0},
    {0,0,D}, {L,0,D},
    {0,H,D}, {L,H,D},
  };

  solid.weights = std::vector<double>(solid.Q.size(),1.0);
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
      writeStoredVar<primary_var_t>("U", [](auto const &x) { return x[1]; });
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
    // for (auto id : sids x.push_back(op(var[id]));

    // add displacement!
    double scale = 1;
    Shape deformed = _shape;
    for (size_t i = 0; i < sids.size(); i++)
    {
      auto const id = sids[i];
      DynamicVectorR const ui = var[id];
      x.push_back(op(ui));
      deformed.Q[i] += scale*convert::to<std::vector<double>>(ui);
    }

      _writer(deformed,x,path+fname+extension);
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

class Stiffness 
{
public:
  using value_t = DynamicMatrixR;
  Stiffness(double E = 70e9, double nu = 0.3)
    : E(E), nu(nu) {}

  ~Stiffness() = default;

  template<typename Point>
  value_t operator()(Point const &p) const
  {
    double const a = E/((1.0+nu)*(1.0-2.0*nu));
    double const b = 1.0-nu;
    double const c = 0.5*(1.0-2.0*nu);

    DynamicMatrixR D(6,6); D.setZero();

    D(0,0) = b;
    D(0,1) = nu;
    D(0,2) = nu;

    D(1,0) = nu;
    D(1,1) = b;
    D(1,2) = nu;

    D(2,0) = nu;
    D(2,1) = nu;
    D(2,2) = b;

    D(3,3) = c;
    D(4,4) = c;
    D(5,5) = c;

    return a*D;
  }

double E,nu;

};

template<size_t NDOF>
struct IdOperator
{
  using value_t = SparseMatrixR;

  template<typename Index>
  value_t const operator()(Index const &p) const
  {
    DynamicVectorR const shape = p.mapper.shape(p.para);
    DynamicMatrixR b(DynamicMatrixR::Zero(NDOF,NDOF*shape.rows()));
    for (int i = 0; i < shape.rows(); i++)
    {
      for (size_t j = 0; j < NDOF; j++)
      {
        b(j,NDOF*i+j) = shape[i];
      }
    }
    return b.sparseView();
  }
};

struct DeformationGradient
{
  using value_t = StaticMatrixR<3,3>;
  using disp_t  = DisplacementVariable<3>;

  DeformationGradient() : _u() {}
  ~DeformationGradient() = default;

  template<typename Index>
  value_t const operator()(Index const &p) const
  {
    auto const grad = p.mapper.grad(p.para);
    auto const eye = value_t::Identity();
    StaticMatrixR<3,3> const Ft = eye + grad*_u(p);
    return Ft.transpose();
  }

  disp_t const _u;
};

class GradOperator 
{
public:
  using value_t = SparseMatrixR;
  GradOperator(): _F() {}
  ~GradOperator() = default;

  template<typename Point>
  value_t operator()(Point const &p) const
  {
    auto const grad = p.mapper.grad(p.para);
    auto const dof = grad.cols();
    auto const F = _F(p);
    DynamicMatrixR b(DynamicMatrixR::Zero(6,3*dof));

    for (int i = 0; i < dof; i++)
    {
      b(0,i*3)   = grad(0,i)*F(0,0);
      b(0,i*3+1) = grad(0,i)*F(1,0);
      b(0,i*3+2) = grad(0,i)*F(2,0);

      b(1,i*3)   = grad(1,i)*F(0,1);
      b(1,i*3+1) = grad(1,i)*F(1,1);
      b(1,i*3+2) = grad(1,i)*F(2,1);

      b(2,i*3)   = grad(2,i)*F(0,2);
      b(2,i*3+1) = grad(2,i)*F(1,2);
      b(2,i*3+2) = grad(2,i)*F(2,2);

      b(3,i*3)   = grad(0,i)*F(0,1)+grad(1,i)*F(0,0);
      b(3,i*3+1) = grad(0,i)*F(1,1)+grad(1,i)*F(1,0);
      b(3,i*3+2) = grad(0,i)*F(2,1)+grad(1,i)*F(2,0);

      b(4,i*3)   = grad(1,i)*F(0,2)+grad(2,i)*F(0,1);
      b(4,i*3+1) = grad(1,i)*F(1,2)+grad(2,i)*F(1,1);
      b(4,i*3+2) = grad(1,i)*F(2,2)+grad(2,i)*F(2,1);

      b(5,i*3)   = grad(0,i)*F(0,2)+grad(2,i)*F(0,0);
      b(5,i*3+1) = grad(0,i)*F(1,2)+grad(2,i)*F(1,0);
      b(5,i*3+2) = grad(0,i)*F(2,2)+grad(2,i)*F(2,0);
    }

    return b.sparseView();
  }

  DeformationGradient const _F;
};

class LinGradOperator 
{
public:
  using value_t = SparseMatrixR;
  LinGradOperator() = default;
  ~LinGradOperator() = default;

  template<typename Point>
  value_t operator()(Point const &p) const
  {
    auto const grad = p.mapper.grad(p.para);
    auto const dof = grad.cols();
    DynamicMatrixR b(DynamicMatrixR::Zero(6,3*dof));

    for (int i = 0; i < dof; i++)
    {
      b(0,i*3)   = grad(0,i);
      b(1,i*3+1) = grad(1,i);
      b(2,i*3+2) = grad(2,i);

      b(3,i*3)   = grad(1,i);
      b(3,i*3+1) = grad(0,i);

      b(4,i*3+1) = grad(2,i);
      b(4,i*3+2) = grad(1,i);

      b(5,i*3)   = grad(2,i);
      b(5,i*3+2) = grad(0,i);
    }

    return b.sparseView();
  }
};

class NLGradOperator 
{
public:
  using value_t = SparseMatrixR;
  NLGradOperator(){}
  ~NLGradOperator() = default;

  template<typename Point>
  value_t operator()(Point const &p) const
  {
    auto const grad = p.mapper.grad(p.para);
    auto const dof = grad.cols();
    DynamicMatrixR b(DynamicMatrixR::Zero(9,3*dof));

    for (int i = 0; i < dof; i++)
    {
      b(0,i*3) = grad(0,i);
      b(1,i*3) = grad(1,i);
      b(2,i*3) = grad(2,i);

      b(3,i*3+1) = grad(0,i);
      b(4,i*3+1) = grad(1,i);
      b(5,i*3+1) = grad(2,i);

      b(6,i*3+2) = grad(0,i);
      b(7,i*3+2) = grad(1,i);
      b(8,i*3+2) = grad(2,i);
    }

    return b.sparseView();
  }
};

struct Strain
{
  using value_t = DynamicVectorR;

  Strain() : _F() {}
  ~Strain() = default;

  template<typename Index>
  value_t const operator()(Index const &p) const
  {
    auto const F = _F(p);
    auto const C = F.transpose()*F;
    auto const eye = StaticMatrixR<3,3>::Identity();
    auto const E = 0.5*(C-eye);

    DynamicVectorR voigt(6);
    voigt[0] = E(0,0);
    voigt[1] = E(1,1);
    voigt[2] = E(2,2);
    voigt[3] = 2.0*E(0,1);
    voigt[4] = 2.0*E(1,2);
    voigt[5] = 2.0*E(0,2);
    return voigt;
  }

  DeformationGradient const _F;
};

struct LinStrain
{
  using value_t = DynamicVectorR;
  using disp_t  = DisplacementVariable<3>;

  LinStrain() : _b(), _u() {}
  ~LinStrain() = default;

  template<typename Index>
  value_t const operator()(Index const &p) const
  {
    DynamicMatrixR const u = _u(p).transpose();
    Eigen::Map<DynamicVectorR const> const uv(u.data(),u.size());
    return _b(p)*uv;
  }

  GradOperator const _b;
  disp_t const _u;
};

struct Stress
{
  using value_t = DynamicVectorR;

  template<typename ...Args>
  Stress(Args const &...args) : _d(args...), _eps() {}
  ~Stress() = default;

  template<typename Index>
  value_t const operator()(Index const &p) const
  {
    return _d(p)*_eps(p);
  }

  Stiffness const _d;
  Strain const _eps;
};

template<typename Stress>
struct ExpandedStress 
{
  using value_t = DynamicMatrixR;

  template<typename ...Args>
  ExpandedStress(Args const &...args) : _stress(args...) {}

  ~ExpandedStress() = default;

  StaticMatrixR<3,3> const fromVoigt(DynamicVectorR const in) const
  {
    StaticMatrixR<3,3> out;
    out(0,0) = in(0);
    out(1,1) = in(1);
    out(2,2) = in(2);
    out(0,1) = in(3);
    out(1,0) = in(3);
    out(1,2) = in(4);
    out(2,1) = in(4);
    out(0,2) = in(5);
    out(2,0) = in(5);
    return out;
  }

  template<typename Index>
  value_t const operator()(Index const &p) const
  {
    auto const stress = fromVoigt(_stress(p));

    DynamicMatrixR out(9,9); out.setZero();
    for (size_t i = 0; i < 3; i++)
    {
      out.block(3*i,3*i,3,3) = stress;
    }
    return out.sparseView();
  }

  Stress const _stress;
};

template<typename T>
class AppliedTraction 
  : public BoundaryConditionBase 
{
public:

private:
  class _f;

public:
  using id_op = IdOperator<3>;
  using linear_form = LinearWeakForm<id_op, _f, DynamicVectorR>;

  template<typename ...Args>
  AppliedTraction(Args ...args)
    : BoundaryConditionBase("AppliedTraction")
    , _load(args...)
    , _load_form(id_op(), _load)
  {}

  virtual ~AppliedTraction() = default;

  WeakFormBase const * getRhs() const override
  { return dynamic_cast<WeakFormBase const *>(&_load_form); }

private:

  class _f
  {
  public:
    using value_t = typename T::value_t;

    template<typename ...Args>
    _f(Args ...args) : _value(args...) {}
    virtual ~_f() {};

    template<typename Index>
    value_t operator()(Index const &p) const
    {
      return _value(p);
    }

  private:
    T _value;
  };

  _f const _load;
  linear_form const _load_form;
};

/*
/////////////////////////////////////////////////////////////////////////////////////
Constraint
/////////////////////////////////////////////////////////////////////////////////////
*/

template<typename ValueType>
class AppliedDisplacement
  : public DirichletConstraint<ValueType, DisplacementVariable<3>>
{
public:
  template<typename ...Args>
  AppliedDisplacement(Args ...args) : DirichletConstraint<ValueType, DisplacementVariable<3>>(args...) {}
  ~AppliedDisplacement() = default;
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
  // NurbsSolid shape; elbow(1.0, 2.0, 3.0, shape);
  double L = 10.0;
  double H = 1.0;
  double D = 1.0;
  NurbsSolid shape; abeam(L,H,D,shape);
  spline_ops::midpointRefinement(3,0,shape);

  // write to file
  // spline_ops::writeToFile(shape,geometryFile,20,2,2);
  // show geometry
  // show(shape,geometryFile);
  {
    double E = 1.2e4; double nu = 0.2;
    // configure the weak forms
    using bdb_form = BdBWeakForm<GradOperator, Stiffness, DynamicMatrixR>;
    using nl_form  = BdBWeakForm<NLGradOperator, ExpandedStress<Stress>, DynamicMatrixR>;
    // add the weak forms
    auto &wfm = WeakFormManager::instance();
    wfm.add<bdb_form>(GradOperator(), Stiffness(E,nu));
    wfm.add<nl_form >(NLGradOperator(), ExpandedStress<Stress>(E,nu));

    // residual
    using residual_form = LinearWeakForm<GradOperator, NegativeValue<Stress>, DynamicVectorR>;
    wfm.add<residual_form>(GradOperator(), NegativeValue<Stress>(E,nu));
  } 
  // add variable storage
  {
    auto mgr = GeometricDofManager::instance();
    mgr.addShape(&shape);
    auto &var = VariableManager::instance();
    auto cpts = mgr.ids.size(); // total of unique cpts
    var.add<primary_var_t>(cpts);
    var.add<OldStoredDisplacementVariable<3>>(cpts);
    var.add<StoredDisplacementVariableIncrement<3>>(cpts);
    var.add<StoredVelocityVariable<3>>(cpts);
    var.add<OldStoredVelocityVariable<3>>(cpts);
  }

  // set up the boundary conditinos
  NurbsSurface fixed; spline_ops::getIsoSurface(0,0,shape,fixed);
  NurbsSurface traction; spline_ops::getIsoSurface(1,1,shape,traction);
  // spline_ops::writeToFile(fixed,geometryFile,2,2);
  // show(fixed,geometryFile);
  // spline_ops::writeToFile(traction,geometryFile,20,2);
  // show(traction,geometryFile);
  {
    using value_t = transpose_type<typename primary_var_t::value_t>::type;
    using fixed_T = AppliedDisplacement<ConstantValue<value_t>>;
    auto &mgr = ConstraintManager::instance();
    mgr.add<fixed_T>(fixed,value_t::Zero());
  }

  {
    using value_t = transpose_type<typename primary_var_t::value_t>::type;
    // using load_T = AppliedTraction<ConstantValue<value_t>>;
    using load_T = AppliedTraction<RampedValue<value_t>>;
    auto &mgr = BoundaryConditionManager::instance();
    mgr.add<load_T>(traction,value_t(0, -32, 0),1.0,value_t::Zero());
  }

  // solve
  {
    using lasystem_t = FeLaSystem<primary_var_t>;
    // using solver_t = LinearSolver<lasystem_t, primary_var_t>;
    // solver_t solver(&shape);
    // using inner_t = LinearSolver<lasystem_t, primary_var_t>;
    using inner_t = NonlinearSolver<lasystem_t, primary_var_t>;
    using solver_t = ContinuationSolver<inner_t>;
    double const tol{1e-6};
    size_t const inner_steps{15};
    size_t const load_increments{12};
    solver_t solver(tol,inner_steps,&shape);
    solver.setStepSize(1.0/load_increments);
    // // add results handling
    using results_t = Results<NurbsSolid,SolidWriter>;
    solver.addResultsHandler<results_t>(shape, SolidWriter(20,2,2), simulation_name);
    solver.run();
  }

  {
    // using results_t = Results<NurbsSolid,SolidWriter>;
    // // show results
    // results_t results(shape, SolidWriter(12,3,24), simulation_name);
    // results.write();
  }
}

int main(int argc, char **argv)
{
  sim();
  return 0;
}
