#pragma once

#include "base/SimulationClock.h"
#include "base/StoredVariableTraits.h"
#include "base/VariableManager.h"
#include "base/VariableUpdateManager.h"

#include "IgaSolverBase.h"
#include "IterativeSolver.h"

#include <algorithm>
#include <cassert>

template<typename LaSystem_t, typename PrimaryVariable>
class NonlinearParabolicInnerSolver
 : public IgaSolverBase, public IterativeSolver
{
public:

  using Variable_t = PrimaryVariable;
  
  template<typename ...Args>
  NonlinearParabolicInnerSolver(double tolerance, size_t max, Args const&...args)
   : _tolerance(tolerance), _maxSteps(max), _lasystem(args...), _error(0.0)
  {}

  ~NonlinearParabolicInnerSolver() = default;

  void preIterate() override
  {
    using u_t = PrimaryVariable;
    using u0_t = typename StoredVariableTraits::old_var<u_t>::type;
    using v0_t = typename StoredVariableTraits::velocity_var<u0_t>::type;

    double alpha1 = _s*_dt;
    double alpha2 = _dt*(1.0-_s);

    // compute fhat
    _lasystem.rhs(_fhat);
    auto &mgr = VariableManager::instance();
    auto u_var = mgr.has<u_t>();
    auto u0_var = mgr.has<u0_t>();
    Eigen::Map<DynamicVectorR> u(&(u_var->data())[0], _dof);
    Eigen::Map<DynamicVectorR> u0(&(u0_var->data())[0], _dof);

    _f = _fhat.head(_dof);
    DynamicVectorR fhat = (_Clast - alpha2*_Klast)*u0 + alpha1*_f + alpha2*_flast - _Clast*u;
    _fhat.head(_dof) = fhat;

    // compute _A
    _lasystem.composed_lhs(_A, 1.0, alpha1);

    // std::cout << "\n***************************" << std::endl;
    // std::cout << "***************************\n" << std::endl;
    // Eigen::IOFormat fmt(4.0);
    // auto f1 = ((1.0/alpha1)*_Clast - alpha3*_Klast)*u0;
    // auto f2 = alpha3*_flast;
    // auto f3 = (1.0/alpha1)*_Clast*u;
    // std::cout << "Sim time = " << SimulationClock::instance().time() << std::endl;
    // std::cout << "Sim step = " << SimulationClock::instance().step() << std::endl;
    // std::cout << "alpha1 = " << alpha1 << std::endl;
    // std::cout << "alpha2 = " << alpha2 << std::endl;
    // std::cout << "alpha3 = " << alpha3 << std::endl;
    // mpm::print(u0_var->data(), "u0 = ");
    // mpm::print(u_var->data(), "u = ");
    // std::cout << "f1 =\n " << f1.transpose().format(fmt) << std::endl;
    // std::cout << "f2 =\n " << f2.transpose().format(fmt) << std::endl;
    // std::cout << "f3 =\n " << f3.transpose().format(fmt) << std::endl;
    // std::cout << "flast =\n " << _flast.transpose().format(fmt) << std::endl;
    // std::cout << "f =\n " << _f.transpose().format(fmt) << std::endl;
    // std::cout << "fhat =\n " << _fhat.transpose().format(fmt) << std::endl;
    // std::cout << "klast =\n " << DynamicMatrixR(_Klast).format(fmt) << std::endl;
    // std::cout << "clast =\n " << DynamicMatrixR((1.0/alpha1)*_Clast).format(fmt) << std::endl;
    // std::cout << "A =\n " << DynamicMatrixR(_A).format(fmt) << std::endl;
    // std::cout << "\n***************************" << std::endl;
    // std::cout << "***************************\n" << std::endl;
  }

  void iterativeSolve() override
  {
    // solve
    Eigen::SparseLU<SparseMatrixR> solver;
    solver.compute(_A);
    if (solver.info() != Eigen::Success)
    {
      std::cout << "Solver decomposition failed!" << std::endl;
      return;
    }

    _dx = solver.solve(_fhat);

    if (solver.info() != Eigen::Success)
    {
      std::cout << "Solver failed!" << std::endl;
      return;
    }
  }

  void postIterate() override
  {
    updateVariableIncrementStorage();
    updateVariableStorage();
    updateVelocityStorage();
    VariableUpdateManager::instance().update();
  }

  void initializeSolution() override
  {
    initializeVariableStorage();
    applyIcs();
    setInitialSystem();
  }

  void solve() override
  {
    _error = 1.0;
    size_t step = 0;
    while (!converged() && step < _maxSteps)
    {
      iterate();
      errorNorm();
      step++;
    }

    if (!converged()) 
    {
      std::cout << "*** INNER SOLVER FAILED TO CONVERGE!!! ***" << std::endl;
      std::cout << "    Error: " << _error << std::endl;
    }
  }

  void finalizeSolution() override
  {
    // no work here
  }

  inline bool converged() const override
  {
    return (_error <= _tolerance);
  }

  void initializeInnerSolve()
  {
    _flast = _f;
    _lasystem.unconstrained_stiffness(_Klast);
    _lasystem.unconstrained_mass(_Clast);
  }

  void setInterpolationFactor(double s) { _s = s; }
  void setTimeStep(double dt) { _dt = dt; }

private:
  void initializeVariableStorage()
  {
    VariableManager::instance().initialize();
    _dof = VariableManager::instance().has<PrimaryVariable>()->size();
  }

  void updateVariableIncrementStorage()
  {
    using inc_t = typename StoredVariableTraits::increment_var<PrimaryVariable>::type;
    auto var = VariableManager::instance().has<inc_t>();
    std::copy(_dx.data(), _dx.data()+var->size(), var->data().begin());
  }

  void updateVariableStorage()
  {
    using inc_t = typename StoredVariableTraits::increment_var<PrimaryVariable>::type;
    using value_t = typename PrimaryVariable::value_t;
    auto var = VariableManager::instance().has<PrimaryVariable>();
    auto dvar = VariableManager::instance().has<inc_t>();
    auto &a= var->data();
    auto &b= dvar->data();
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::plus<value_t>());
  }

  void updateVelocityStorage()
  {
    // store velocity variable
    using u_t = PrimaryVariable;
    using v_t =  typename StoredVariableTraits::velocity_var<u_t>::type;
    using u0_t = typename StoredVariableTraits::old_var<u_t>::type;
    using v0_t = typename StoredVariableTraits::old_var<v_t>::type;

    auto &mgr = VariableManager::instance();

    auto u_var = mgr.has<u_t>(); 
    auto v_var = mgr.has<v_t>(); 
    auto u0_var = mgr.has<u0_t>(); 
    auto v0_var = mgr.has<v0_t>(); 

    auto dof = u_var->size();    
    Eigen::Map<DynamicVectorR>  u(&(u_var->data())[0], dof);
    Eigen::Map<DynamicVectorR> u0(&(u0_var->data())[0], dof);
    Eigen::Map<DynamicVectorR> v0(&(v0_var->data())[0], dof);

    double alpha1 = 1.0/(_s*_dt);
    double alpha2 = (1.0 - _s)/_s;
    DynamicVectorR v = alpha1*(u - u0) - alpha2*v0;
    std::copy(v.data(), v.data()+dof,  v_var->data().begin());
  }
  
  void applyIcs()
  {
    std::vector<Triplet> initialbc;
    auto const &constrained = _lasystem.constrained_system();
    constrained.sparse_rhs(initialbc);

    using u_t = PrimaryVariable;
    auto primary = VariableManager::instance().has<u_t>();
    auto &var = primary->data();
    for_each(initialbc.begin(), initialbc.end(),
      [&](Triplet const &t) { var[t.row()] = t.value(); }
      );
  }

  void setInitialSystem()
  {
    _lasystem.unconstrained_stiffness(_Klast);
    _lasystem.unconstrained_mass(_Clast);
    _flast.resize(_dof);
    _f.resize(_dof);
    std::fill(_flast.data(), _flast.data()+_dof, 0.0);
    std::fill(_f.data(), _f.data()+_dof, 0.0);
  }

  inline void errorNorm()
  {
    _error = _dx.head(_dof).norm()/_dof;
  }

  double _tolerance;
  size_t _maxSteps;
  LaSystem_t const _lasystem;
  double _error;
  double _s;
  double _dt;
  size_t _dof;
  SparseMatrixR _A, _Klast, _Clast;
  DynamicVectorR _dx, _f, _flast, _fhat;
};