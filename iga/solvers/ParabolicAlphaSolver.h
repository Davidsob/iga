#pragma once

#include "base/StoredVariableTraits.h"
#include "base/VariableManager.h"

#include "LinearParabolicInnerSolver.h"
#include "NonlinearParabolicInnerSolver.h"
#include "TimeDependentSolver.h"

#include <map>
#include <algorithm>

template<typename InnerSolver_t>
class ParabolicAlphaSolver_T
  : public TimeDependentSolver<InnerSolver_t>
{
public:

  using base_t = TimeDependentSolver<InnerSolver_t>;

  template<typename ...Args>
  ParabolicAlphaSolver_T(Args const&...args)
    : base_t(args...)
    , _t0(0), _tf(0), _dt(1), _step(0)
  {}

  ~ParabolicAlphaSolver_T() = default;

  void setStartTime(double t0) { _t0 = t0; }
  void setFinalTime(double tf) { _tf = tf; }
  void setStepSize(double dt) {
    _dt = dt;
    this->_innerSolver.setTimeStep(dt);
  }
  void setInterpolationFactor(double s) {
    this->_innerSolver.setInterpolationFactor(s);
  }

  double startTime() const override
  {
    return _t0;
  }

  double finishTime() const override
  {
    return _tf;
  }

  void initializeInnerSolve()
  {
    this->_innerSolver.initializeInnerSolve();
    // store old vars
    using u_t = typename InnerSolver_t::Variable_t;
    using v_t =  typename StoredVariableTraits::velocity_var<u_t>::type;
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
    SimulationClock::instance() += _dt;
    _step ++;
  }

  void initializeSolution() override final
  {
    base_t::initializeSolution();
    this->_innerSolver.initializeSolution();
  }

  void finalizeSolution() override final
  {
    //
  }

private:
  double _t0;
  double _tf;
  double _dt;
  size_t _step;
};

template<typename LaSystem_t, typename PrimaryVariable>
using NLParabolicAlphaSolver = ParabolicAlphaSolver_T<NonlinearParabolicInnerSolver<LaSystem_t, PrimaryVariable>>;

template<typename LaSystem_t, typename PrimaryVariable>
using ParabolicAlphaSolver = ParabolicAlphaSolver_T<LinearParabolicInnerSolver<LaSystem_t, PrimaryVariable>>;