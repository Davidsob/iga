#pragma once

#include "base/SimulationClock.h"
#include "base/SimulationResults.h"
#include "base/StoredVariableTraits.h"

#include "SolverBase.h"

#include <memory>

template<typename InnerSolver>
class TimeDependentSolver
  : public SolverBase
{
public:

  virtual ~TimeDependentSolver() = default;
 
  virtual double startTime() const = 0; 
  virtual double finishTime() const = 0; 
  double currentTime() const
  { return SimulationClock::instance().time(); }

  virtual void initializeInnerSolve() = 0; 
  virtual void finalizeInnerSolve() = 0; 
  virtual void innerSolve()
  { _innerSolver.solve(); }

  virtual void initializeSolution() override
  {
    SimulationClock::instance().set(startTime());
  }

  void solve() override
  {
    while ( currentTime() <= finishTime() && converged())
    {
      initializeInnerSolve();
      innerSolve();
      for (auto const &hndlr : _results) hndlr->write();
      finalizeInnerSolve();
    }
  }

  virtual void finalizeSolution() override {};

  template<typename T, typename ...Args>
  void addResultsHandler(Args const&...args)
  {
    _results.push_back(std::move(std::unique_ptr<SimulationResults>(new T(args...))));
  }

  bool takeStep()
  {
    initializeInnerSolve();
    innerSolve();
    for (auto const &hndlr : _results) hndlr->write();
    finalizeInnerSolve();
    return converged();
  }

  bool converged() const { return _innerSolver.converged(); }

protected:

  template<typename ...Args>
  TimeDependentSolver(Args const&...args)
    : _innerSolver(args...)
  {}

  InnerSolver _innerSolver;

private:
  using results_t = std::unique_ptr<SimulationResults>;
   std::vector<results_t> _results;
};