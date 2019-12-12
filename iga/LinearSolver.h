#pragma once

#include "base/SolverBase.h"
#include "base/SimulationClock.h"
#include "base/VariableManager.h"
#include "base/VariableTags.h"
#include "base/VariableUpdateManager.h"

#include <Eigen/SparseLU>

template<typename LaSystem_t, typename PrimaryVariable>
class LinearSolver
  : public SolverBase
{
public:
  template<typename ...Args>
  LinearSolver(Args const &...args)
    : _lasystem(args...)
  {}

  ~LinearSolver() = default;

  void initializeSolution() override
  {
    SimulationClock::instance().set(0.0);
    initializeVariableStorage();
  }

  void solve() override
  {
    std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
    std::cout << "Compile..." << std::endl;
    SparseMatrixR A;
    DynamicVectorR b;
    _lasystem.discretize(A, b);

    // solve
    std::cout << "Begin linear solve..." << std::endl;
    Eigen::SparseLU<SparseMatrixR> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success)
    {
      std::cout << "Solver decomposition failed!" << std::endl;
      return;
    }

    _x.resize(b.size());
    _x = solver.solve(b);
    if (solver.info() != Eigen::Success)
    {
      std::cout << "Solver failed!" << std::endl;
      return;
    }

    std::cout << "Completed linear solve!" << std::endl;
    std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
  }

  void finalizeSolution() override
  {
    updateVariableStorage();
    VariableUpdateManager::instance().update();
  }

private:

  void initializeVariableStorage()
  {
    VariableManager::instance().initialize();
  }

  void updateVariableStorage()
  {
    auto var = VariableManager::instance().has<PrimaryVariable>();
    std::copy(_x.data(), _x.data()+var->size(), var->data().begin());
  }

  LaSystem_t const _lasystem;
  DynamicVectorR _x;
};
