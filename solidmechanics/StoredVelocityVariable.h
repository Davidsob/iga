#pragma once

#include "base/StoredVariable.h"
#include "base/StoredVariableOperator.h"
#include "base/VariableTags.h"

#include "utils/MatrixTypes.h"

template<size_t dim>
class StoredVelocityVariable
  : public StoredVariable<StaticRowVectorR<dim>>
{
public:
  using value_t = StaticRowVectorR<dim>;
  static constexpr size_t const &ndof = dim;

  explicit StoredVelocityVariable(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "StoredVelocityVariable")
  {}

  ~StoredVelocityVariable() = default;
};

template<size_t dim>
class StoredVelocityVariableIncrement
  : public StoredVariable<StaticRowVectorR<dim>>
{
public:
  using value_t = StaticRowVectorR<dim>;
  static constexpr size_t const &ndof = dim;

  explicit StoredVelocityVariableIncrement(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "StoredVelocityVariableIncrement")
  {}

  ~StoredVelocityVariableIncrement() = default;
};

template<size_t dim>
class OldStoredVelocityVariable
  : public StoredVariable<StaticRowVectorR<dim>>
{
public:
  using value_t = StaticRowVectorR<dim>;
  static constexpr size_t const &ndof = dim;

  explicit OldStoredVelocityVariable(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "OldStoredVelocityVariable")
  {}

  ~OldStoredVelocityVariable() = default;
};

template<size_t dim>
class VelocityVariable
  : public StoredVariableOperator<StoredVelocityVariable<dim>>
{
public:
  static constexpr size_t const &ndof = dim;

  VelocityVariable() = default;
  ~VelocityVariable() = default;
};

template<size_t dim>
class VelocityVariableIncrement
  : public StoredVariableOperator<StoredVelocityVariableIncrement<dim>>
{
public:
  static constexpr size_t const &ndof = dim;

  VelocityVariableIncrement() = default;
  ~VelocityVariableIncrement() = default;
};

template<size_t dim>
class OldVelocityVariable
  : public StoredVariableOperator<OldStoredVelocityVariable<dim>>
{
public:
  static constexpr size_t const &ndof = dim;

  OldVelocityVariable() = default;
  ~OldVelocityVariable() = default;
};
