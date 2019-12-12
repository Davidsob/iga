#pragma once

#include "StoredVariable.h"
#include "StoredVariableOperator.h"
#include "VariableTags.h"

class StoredTemperatureRateVariable
  : public StoredVariable<double>
{
public:
  using value_t = double;
  static constexpr size_t const &ndof = 1;

  explicit StoredTemperatureRateVariable(size_t i)  
    : StoredVariable<double>(i, 0, "StoredTemperatureRateVariable")
  {}

  ~StoredTemperatureRateVariable() = default;
};

class StoredTemperatureRateVariableIncrement
  : public StoredVariable<double>
{
public:
  using value_t = double;
  static constexpr size_t const &ndof = 1;

  explicit StoredTemperatureRateVariableIncrement(size_t i)  
    : StoredVariable<double>(i, 0, "StoredTemperatureRateVariableIncrement")
  {}

  ~StoredTemperatureRateVariableIncrement() = default;
};

class OldStoredTemperatureRateVariable
  : public StoredVariable<double>
{
public:
  using value_t = double;
  static constexpr size_t const &ndof = 1;

  explicit OldStoredTemperatureRateVariable(size_t i)  
    : StoredVariable<double>(i, 0, "OldStoredTemperatureRateVariable")
  {}

  ~OldStoredTemperatureRateVariable() = default;
};

class TemperatureRateVariable
  : public StoredVariableOperator<StoredTemperatureRateVariable>
{
public:
  static constexpr size_t const &ndof = 1;
  
  TemperatureRateVariable() = default;
  ~TemperatureRateVariable() = default;
};

class TemperatureRateVariableIncrement
  : public StoredVariableOperator<StoredTemperatureRateVariableIncrement>
{
public:
  static constexpr size_t const &ndof = 1;
  
  TemperatureRateVariableIncrement() = default;
  ~TemperatureRateVariableIncrement() = default;
};

class OldTemperatureRateVariable
  : public StoredVariableOperator<OldStoredTemperatureRateVariable>
{
public:
  static constexpr size_t const &ndof = 1;
  
  OldTemperatureRateVariable() = default;
  ~OldTemperatureRateVariable() = default;
};