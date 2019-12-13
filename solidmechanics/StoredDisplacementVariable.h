#pragma once

#include "base/StoredVariable.h"
#include "base/StoredVariableOperator.h"
#include "base/VariableTags.h"

#include "utils/MatrixTypes.h"

template<size_t dim>
class StoredDisplacementVariable
  : public StoredVariable<StaticRowVectorR<dim>>
{
public:
  using value_t = StaticRowVectorR<dim>;
  static constexpr size_t const &ndof = dim;

  explicit StoredDisplacementVariable(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "StoredDisplacementVariable")
  {

  }

  ~StoredDisplacementVariable() = default;
};

template<size_t dim>
class StoredDisplacementVariableIncrement
  : public StoredVariable<StaticRowVectorR<dim>>
{
public:
  using value_t = StaticRowVectorR<dim>;
  static constexpr size_t const &ndof = dim;

  explicit StoredDisplacementVariableIncrement(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "StoredDisplacementVariableIncrement")
  {}

  ~StoredDisplacementVariableIncrement() = default;
};

template<size_t dim>
class OldStoredDisplacementVariable
  : public StoredVariable<StaticRowVectorR<dim>>
{
public:
  using value_t = StaticRowVectorR<dim>;
  static constexpr size_t const &ndof = dim;

  explicit OldStoredDisplacementVariable(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "OldStoredDisplacementVariable")
  {}

  ~OldStoredDisplacementVariable() = default;
};

template<size_t dim>
class DisplacementVariable
  : public StoredVariableOperator<StoredDisplacementVariable<dim>>
{
public:
  static constexpr size_t const &ndof = dim;

  DisplacementVariable() = default;
  ~DisplacementVariable() = default;
};

template<size_t dim>
class DisplacementVariableIncrement
  : public StoredVariableOperator<StoredDisplacementVariableIncrement<dim>>
{
public:
  static constexpr size_t const &ndof = dim;

  DisplacementVariableIncrement() = default;
  ~DisplacementVariableIncrement() = default;
};

template<size_t dim>
class OldDisplacementVariable
  : public StoredVariableOperator<OldStoredDisplacementVariable<dim>>
{
public:
  static constexpr size_t const &ndof = dim;

  OldDisplacementVariable() = default;
  ~OldDisplacementVariable() = default;
};
