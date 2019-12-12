#pragma once

#include "base/SimulationClock.h"
#include "base/VariableManager.h"

#include "utils/MatrixTypes.h"

#include "TimeDependentSolver.h"

template<typename InnerSolver>
class ContinuationSolver
  : public TimeDependentSolver<InnerSolver>
{
public:
  template<typename ...Args>
  ContinuationSolver(Args const&...args)
    : TimeDependentSolver<InnerSolver>(args...)
    , _dt(1.0)
  {}

  ~ContinuationSolver() = default;

  double startTime() const override
  {
    return 0;
  }

  double finishTime() const override
  {
    return 1.0;
  }

  void setStepSize(double dt) { _dt = dt; }

  void initializeInnerSolve()
  {
    // update simulation time
    using u_t  = typename InnerSolver::Variable_t;
    using v_t  = typename StoredVariableTraits::velocity_var<u_t>::type;
    using u0_t = typename StoredVariableTraits::old_var<u_t>::type;
    using v0_t = typename StoredVariableTraits::old_var<v_t>::type;

    auto &mgr = VariableManager::instance();

    auto u_var = mgr.has<u_t>(); 
    auto v_var = mgr.has<v_t>(); 
    auto u0_var = mgr.has<u0_t>(); 
    auto v0_var = mgr.has<v0_t>(); 
    u0_var->data() = u_var->data();
    v0_var->data() = v_var->data();
  }

  void finalizeInnerSolve()
  {
    // iterate simulation time
    // std::cout << "*** Solution for t [s] = " << SimulationClock::instance().time() << std::endl;
    // using var_t = typename InnerSolver::Variable_t;
    // auto var = VariableManager::instance().has<var_t>();
    // auto &x = var->data();
    // std::cout << "   " << Eigen::Map<DynamicVectorR>(&x[0], x.size()) << std::endl;
    SimulationClock::instance() += _dt;
  }

  void initializeSolution() override final
  {
    TimeDependentSolver<InnerSolver>::initializeSolution();
    this->_innerSolver.initializeSolution();
  }

  void finalizeSolution() override final
  {
    //
  }

private:
  double _dt;
};