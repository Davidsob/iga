#pragma once

#include "base/StoredVariable.h"
#include "base/StoredVariableOperator.h"
#include "base/VariableTags.h"

class StoredTemperatureVariable
  : public StoredVariable<double>
{
public:
  using value_t = double;
  static constexpr size_t const &ndof = 1;

  explicit StoredTemperatureVariable(size_t i)  
    : StoredVariable<double>(i, 0, "StoredTemperatureVariable")
  {}

  ~StoredTemperatureVariable() = default;
};

class StoredTemperatureVariableIncrement
  : public StoredVariable<double>
{
public:
  using value_t = double;
  static constexpr size_t const &ndof = 1;

  explicit StoredTemperatureVariableIncrement(size_t i)  
    : StoredVariable<double>(i, 0, "StoredTemperatureVariableIncrement")
  {}

  ~StoredTemperatureVariableIncrement() = default;
};

class OldStoredTemperatureVariable
  : public StoredVariable<double>
{
public:
  using value_t = double;
  static constexpr size_t const &ndof = 1;

  explicit OldStoredTemperatureVariable(size_t i)  
    : StoredVariable<double>(i, 0, "OldStoredTemperatureVariable")
  {}

  ~OldStoredTemperatureVariable() = default;
};

class TemperatureVariable
  : public StoredVariableOperator<StoredTemperatureVariable>
{
public:
  static constexpr size_t const &ndof = 1;

  TemperatureVariable() = default;
  ~TemperatureVariable() = default;
};

class TemperatureVariableIncrement
  : public StoredVariableOperator<StoredTemperatureVariableIncrement>
{
public:
  static constexpr size_t const &ndof = 1;

  TemperatureVariableIncrement() = default;
  ~TemperatureVariableIncrement() = default;
};

class OldTemperatureVariable
  : public StoredVariableOperator<OldStoredTemperatureVariable>
{
public:
  static constexpr size_t const &ndof = 1;

  OldTemperatureVariable() = default;
  ~OldTemperatureVariable() = default;
};
