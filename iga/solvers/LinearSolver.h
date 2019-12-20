#pragma once

#include "base/SimulationClock.h"
#include "base/VariableManager.h"
#include "base/VariableTags.h"
#include "base/VariableUpdateManager.h"

#include "splines/utils/Converters.h"

#include <Eigen/SparseLU>

#include "IgaSolverBase.h"

template<typename LaSystem_t, typename PrimaryVariable>
class LinearSolver
  : public IgaSolverBase
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
    std::printf("Lhs(%lu,%lu)\n",A.rows(), A.cols());
    std::printf("Rhs(%lu)\n",b.size());
    // IO::spy(DynamicMatrixR(A));
    Eigen::SparseLU<SparseMatrixR> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success)
    {
      std::cout << "Solver decomposition failed!" << std::endl;
      return;
    }

    _x.resize(b.size());
    std::cout << "Sum f = " << b.sum() << std::endl;
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
    using to_t = std::vector<typename PrimaryVariable::value_t>;
    auto const ndof = PrimaryVariable::ndof;
    auto const dof  = _x.size();
    DynamicMatrixR const xreshaped = Eigen::Map<DynamicMatrixR const>(_x.data(),ndof,dof/ndof).transpose();

    auto const x = convert::to<to_t>(xreshaped);
    auto var = VariableManager::instance().has<PrimaryVariable>();
    std::copy(x.begin(), x.begin()+var->size(), var->data().begin());
  }

  LaSystem_t const _lasystem;
  DynamicVectorR _x;
};
