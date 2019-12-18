#pragma once

#include "base/StoredVariableTraits.h"
#include "base/VariableUpdateManager.h"

#include "splines/utils/Converters.h"

#include "IgaSolverBase.h"
#include "IterativeSolver.h"

#include <algorithm>

template<typename LaSystem_t, typename PrimaryVariable>
class NonlinearSolver
 : public IgaSolverBase, public IterativeSolver
{
public:

  using Variable_t = PrimaryVariable;
  
  template<typename ...Args>
  NonlinearSolver(double tolerance, size_t max, Args const&...args)
   : _tolerance(tolerance), _maxSteps(max), _lasystem(args...)
  {}

  ~NonlinearSolver() = default;

  void preIterate() override
  {

  }

  void iterativeSolve() override
  {
    SparseMatrixR A;
    DynamicVectorR b;
    _lasystem.discretize(A, b);

    // solve
    Eigen::SparseLU<SparseMatrixR> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success)
    {
      std::cout << "Solver decomposition failed!" << std::endl;
      return;
    }

    _dx.resize(b.size());
    _dx = solver.solve(b);
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
    VariableUpdateManager::instance().update();
  }

  void initializeSolution() override
  {
    initializeVariableStorage();
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
    } else {
      std::cout << "*** INNER SOLVER CONVERGED ***" << std::endl;
      std::cout << "    Steps: " << step   << std::endl;
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

private:
  void initializeVariableStorage()
  {
    VariableManager::instance().initialize();
    _dof = VariableManager::instance().has<PrimaryVariable>()->size();
  }

  void updateVariableIncrementStorage()
  {
    using inc_t = typename StoredVariableTraits::increment_var<PrimaryVariable>::type;
    using to_t = std::vector<typename inc_t::value_t>;
    // std::copy(_dx.data(), _dx.data()+dvar->size(), dvar->data().begin());

    auto const ndof = PrimaryVariable::ndof;
    auto const dof  = _dx.size();
    DynamicMatrixR const dxreshaped = Eigen::Map<DynamicMatrixR const>(_dx.data(),ndof,dof/ndof).transpose();

    auto const x = convert::to<to_t>(dxreshaped);
    auto dvar = VariableManager::instance().has<inc_t>();
    std::copy(x.begin(), x.begin()+dvar->size(), dvar->data().begin());
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

  inline void errorNorm()
  {
    _error = _dx.head(_dof).norm()/_dof;
  }

  double _tolerance;
  size_t _maxSteps;
  LaSystem_t const _lasystem;
  double _error;
  double _dof;
  DynamicVectorR _dx;
};