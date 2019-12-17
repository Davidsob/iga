#pragma once

#include "base/StoredVariable.h"
#include "base/StoredVariableOperator.h"
#include "base/VariableTags.h"

#include "utils/MatrixTypes.h"


class StoredSixDofVelocityVariable
  : public StoredVariable<StaticRowVectorR<6>>
{
public:
  using value_t = StaticRowVectorR<6>;
  static constexpr size_t const &ndof = 6;

  explicit StoredSixDofVelocityVariable(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "StoredSixDofVelocityVariable")
  {}

  ~StoredSixDofVelocityVariable() = default;
};


class StoredSixDofVelocityVariableIncrement
  : public StoredVariable<StaticRowVectorR<6>>
{
public:
  using value_t = StaticRowVectorR<6>;
  static constexpr size_t const &ndof = 6;

  explicit StoredSixDofVelocityVariableIncrement(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "StoredSixDofVelocityVariableIncrement")
  {}

  ~StoredSixDofVelocityVariableIncrement() = default;
};


class OldStoredSixDofVelocityVariable
  : public StoredVariable<StaticRowVectorR<6>>
{
public:
  using value_t = StaticRowVectorR<6>;
  static constexpr size_t const &ndof = 6;

  explicit OldStoredSixDofVelocityVariable(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "OldStoredSixDofVelocityVariable")
  {}

  ~OldStoredSixDofVelocityVariable() = default;
};


class SixDofVelocityVariable
  : public StoredVariableOperator<StoredSixDofVelocityVariable>
{
public:
  static constexpr size_t const &ndof = 6;

  SixDofVelocityVariable() = default;
  ~SixDofVelocityVariable() = default;
};


class SixDofVelocityVariableIncrement
  : public StoredVariableOperator<StoredSixDofVelocityVariableIncrement>
{
public:
  static constexpr size_t const &ndof = 6;

  SixDofVelocityVariableIncrement() = default;
  ~SixDofVelocityVariableIncrement() = default;
};


class OldSixDofVelocityVariable
  : public StoredVariableOperator<OldStoredSixDofVelocityVariable>
{
public:
  static constexpr size_t const &ndof = 6;

  OldSixDofVelocityVariable() = default;
  ~OldSixDofVelocityVariable() = default;
};
