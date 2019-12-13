#pragma once

#include "base/StoredVariable.h"
#include "base/StoredVariableOperator.h"
#include "base/VariableTags.h"

#include "utils/MatrixTypes.h"

template<size_t dim>
class StoredAccelerationVariable
  : public StoredVariable<StaticRowVectorR<dim>>
{
public:
  using value_t = StaticRowVectorR<dim>;
  static constexpr size_t const &ndof = dim;

  explicit StoredAccelerationVariable(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "StoredAccelerationVariable")
  {}

  ~StoredAccelerationVariable() = default;
};

template<size_t dim>
class StoredAccelerationVariableIncrement
  : public StoredVariable<StaticRowVectorR<dim>>
{
public:
  using value_t = StaticRowVectorR<dim>;
  static constexpr size_t const &ndof = dim;

  explicit StoredAccelerationVariableIncrement(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "StoredAccelerationVariableIncrement")
  {}

  ~StoredAccelerationVariableIncrement() = default;
};

template<size_t dim>
class OldStoredAccelerationVariable
  : public StoredVariable<StaticRowVectorR<dim>>
{
public:
  using value_t = StaticRowVectorR<dim>;
  static constexpr size_t const &ndof = dim;

  explicit OldStoredAccelerationVariable(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "OldStoredAccelerationVariable")
  {}

  ~OldStoredAccelerationVariable() = default;
};

template<size_t dim>
class AccelerationVariable
  : public StoredVariableOperator<StoredAccelerationVariable<dim>>
{
public:
  static constexpr size_t const &ndof = dim;

  AccelerationVariable() = default;
  ~AccelerationVariable() = default;
};

template<size_t dim>
class AccelerationVariableIncrement
  : public StoredVariableOperator<StoredAccelerationVariableIncrement<dim>>
{
public:
  static constexpr size_t const &ndof = dim;

  AccelerationVariableIncrement() = default;
  ~AccelerationVariableIncrement() = default;
};

template<size_t dim>
class OldAccelerationVariable
  : public StoredVariableOperator<OldStoredAccelerationVariable<dim>>
{
public:
  static constexpr size_t const &ndof = dim;

  OldAccelerationVariable() = default;
  ~OldAccelerationVariable() = default;
};
