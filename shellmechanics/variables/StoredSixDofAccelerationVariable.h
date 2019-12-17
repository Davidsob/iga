#pragma once

#include "base/StoredVariable.h"
#include "base/StoredVariableOperator.h"
#include "base/VariableTags.h"

#include "utils/MatrixTypes.h"


class StoredSixDofAccelerationVariable
  : public StoredVariable<StaticRowVectorR<6>>
{
public:
  using value_t = StaticRowVectorR<6>;
  static constexpr size_t const &ndof = 6;

  explicit StoredSixDofAccelerationVariable(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "StoredSixDofAccelerationVariable")
  {}

  ~StoredSixDofAccelerationVariable() = default;
};


class StoredSixDofAccelerationVariableIncrement
  : public StoredVariable<StaticRowVectorR<6>>
{
public:
  using value_t = StaticRowVectorR<6>;
  static constexpr size_t const &ndof = 6;

  explicit StoredSixDofAccelerationVariableIncrement(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "StoredSixDofAccelerationVariableIncrement")
  {}

  ~StoredSixDofAccelerationVariableIncrement() = default;
};


class OldStoredSixDofAccelerationVariable
  : public StoredVariable<StaticRowVectorR<6>>
{
public:
  using value_t = StaticRowVectorR<6>;
  static constexpr size_t const &ndof = 6;

  explicit OldStoredSixDofAccelerationVariable(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "OldStoredSixDofAccelerationVariable")
  {}

  ~OldStoredSixDofAccelerationVariable() = default;
};


class SixDofAccelerationVariable
  : public StoredVariableOperator<StoredSixDofAccelerationVariable>
{
public:
  static constexpr size_t const &ndof = 6;

  SixDofAccelerationVariable() = default;
  ~SixDofAccelerationVariable() = default;
};


class SixDofAccelerationVariableIncrement
  : public StoredVariableOperator<StoredSixDofAccelerationVariableIncrement>
{
public:
  static constexpr size_t const &ndof = 6;

  SixDofAccelerationVariableIncrement() = default;
  ~SixDofAccelerationVariableIncrement() = default;
};


class OldSixDofAccelerationVariable
  : public StoredVariableOperator<OldStoredSixDofAccelerationVariable>
{
public:
  static constexpr size_t const &ndof = 6;

  OldSixDofAccelerationVariable() = default;
  ~OldSixDofAccelerationVariable() = default;
};
