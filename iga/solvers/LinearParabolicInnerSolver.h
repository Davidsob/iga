#pragma once

#include "base/SimulationClock.h"
#include "base/StoredVariableTraits.h"
#include "base/VariableManager.h"
#include "base/VariableUpdateManager.h"

#include "SolverBase.h"
#include "IterativeSolver.h"

#include <algorithm>
#include <cassert>

template<typename LaSystem_t, typename PrimaryVariable>
class LinearParabolicInnerSolver
 : public SolverBase
{
public:

  using Variable_t = PrimaryVariable;
  
  template<typename ...Args>
  LinearParabolicInnerSolver(Args const&...args)
   : _lasystem(args...)
  {}

  ~LinearParabolicInnerSolver() = default;

  void compileSystem()
  {
    using u0_t = typename StoredVariableTraits::old_var<PrimaryVariable>::type;
    using v0_t = typename StoredVariableTraits::velocity_var<u0_t>::type;

    double alpha1 = _s*_dt;
    double alpha2 = _dt*(1.0-_s);

    // compute fhat
    _lasystem.rhs(_fhat);
    auto &mgr = VariableManager::instance();
    auto u0_var = mgr.has<u0_t>();
    Eigen::Map<DynamicVectorR> u0(&(u0_var->data())[0], _dof);

    _f = _fhat.head(_dof);
    DynamicVectorR fhat = (_Clast - alpha2*_Klast)*u0 + alpha1*_f + alpha2*_flast;
    _fhat.head(_dof) = fhat;

    // compute _A
    _lasystem.composed_lhs(_A, 1.0, alpha1);
  }

  void solveSystem()
  {
    // solve
    Eigen::SparseLU<SparseMatrixR> solver;
    solver.compute(_A);
    if (solver.info() != Eigen::Success)
    {
      std::cout << "Solver decomposition failed!" << std::endl;
      return;
    }

    _x = solver.solve(_fhat);

    if (solver.info() != Eigen::Success)
    {
      std::cout << "Solver failed!" << std::endl;
      return;
    }
  }

  void updateSolutionVariables()
  {
    updateVariableStorage();
    updateVelocityStorage();
    VariableUpdateManager::instance().update();
    _flast = _f;
    _lasystem.unconstrained_stiffness(_Klast);
    _lasystem.unconstrained_mass(_Clast);
  }

  void initializeSolution() override
  {
    initializeVariableStorage();
    applyIcs();
    setInitialSystem();
  }

  void solve() override
  {
    compileSystem();
    solveSystem();
    updateSolutionVariables();
  }

  void finalizeSolution() override
  {
    // no work here
  }

  inline bool converged() const
  {
    return true; // linear solver always converged
  }

  void setInterpolationFactor(double s) { _s = s; }
  void setTimeStep(double dt) { _dt = dt; }

private:
  void initializeVariableStorage()
  {
    VariableManager::instance().initialize();
    _dof = VariableManager::instance().has<PrimaryVariable>()->size();
  }

  void updateVariableStorage()
  {
    using inc_t = typename StoredVariableTraits::increment_var<PrimaryVariable>::type;
    using value_t = typename PrimaryVariable::value_t;
    auto var = VariableManager::instance().has<PrimaryVariable>();
    std::copy(_x.data(), _x.data()+var->size(), var->data().begin());
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
    auto primary = VariableManager::instance().has<PrimaryVariable>();
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
    std::fill(_flast.data(), _flast.data()+_dof, 0.0);
  }

  LaSystem_t const _lasystem;
  double _s;
  double _dt;
  size_t _dof;
  SparseMatrixR _A, _Klast, _Clast;
  DynamicVectorR _x, _f, _flast, _fhat;
};