#pragma once

#include "base/StoredVariable.h"
#include "base/StoredVariableOperator.h"
#include "base/VariableTags.h"

#include "utils/MatrixTypes.h"

class StoredSixDofDisplacementVariable
  : public StoredVariable<StaticRowVectorR<6>>
{
public:
  using value_t = StaticRowVectorR<6>;
  static constexpr size_t const &ndof = 6;

  explicit StoredSixDofDisplacementVariable(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "StoredSixDofDisplacementVariable")
  {

  }

  ~StoredSixDofDisplacementVariable() = default;
};

class StoredSixDofDisplacementVariableIncrement
  : public StoredVariable<StaticRowVectorR<6>>
{
public:
  using value_t = StaticRowVectorR<6>;
  static constexpr size_t const &ndof = 6;

  explicit StoredSixDofDisplacementVariableIncrement(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "StoredSixDofDisplacementVariableIncrement")
  {}

  ~StoredSixDofDisplacementVariableIncrement() = default;
};

class OldStoredSixDofDisplacementVariable
  : public StoredVariable<StaticRowVectorR<6>>
{
public:
  using value_t = StaticRowVectorR<6>;
  static constexpr size_t const &ndof = 6;

  explicit OldStoredSixDofDisplacementVariable(size_t i)  
    : StoredVariable<value_t>(i, value_t::Zero(), "OldStoredSixDofDisplacementVariable")
  {}

  ~OldStoredSixDofDisplacementVariable() = default;
};

class SixDofDisplacementVariable
  : public StoredVariableOperator<StoredSixDofDisplacementVariable>
{
public:
  static constexpr size_t const &ndof = 6;

  SixDofDisplacementVariable() = default;
  ~SixDofDisplacementVariable() = default;
};

class SixDofDisplacementVariableIncrement
  : public StoredVariableOperator<StoredSixDofDisplacementVariableIncrement>
{
public:
  static constexpr size_t const &ndof = 6;

  SixDofDisplacementVariableIncrement() = default;
  ~SixDofDisplacementVariableIncrement() = default;
};

class OldSixDofDisplacementVariable
  : public StoredVariableOperator<OldStoredSixDofDisplacementVariable>
{
public:
  static constexpr size_t const &ndof = 6;

  OldSixDofDisplacementVariable() = default;
  ~OldSixDofDisplacementVariable() = default;
};
