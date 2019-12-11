#pragma once

#include "StoredVoltageVariable.h"
#include "StoredElectricVelocityVariable.h"

#include "StoredTemperatureVariable.h"
#include "StoredTemperatureVelocityVariable.h"

namespace StoredVariableTraits
{

  template<typename Var> struct increment_var {};

  // template<>
  // struct increment_var<StoredVoltageVariable>
  // {
  //   using type = StoredVoltageVariableIncrement;
  // };

  // template<>
  // struct increment_var<StoredElectricVelocityVariable>
  // {
  //   using type = StoredElectricVelocityVariableIncrement;
  // };

  // template<>
  // struct increment_var<StoredTemperatureVariable>
  // {
  //   using type = StoredTemperatureVariableIncrement;
  // };

  // template<>
  // struct increment_var<StoredTemperatureVelocityVariable>
  // {
  //   using type = StoredTemperatureVelocityVariableIncrement;
  // };

  template<typename Var> struct old_var {};

  // template<>
  // struct old_var<StoredVoltageVariable>
  // {
  //   using type = OldStoredVoltageVariable;
  // };

  // template<>
  // struct old_var<StoredElectricVelocityVariable>
  // {
  //   using type = OldStoredElectricVelocityVariable;
  // };

  // template<>
  // struct old_var<StoredTemperatureVariable>
  // {
  //   using type = OldStoredTemperatureVariable;
  // };

  // template<>
  // struct old_var<StoredTemperatureVelocityVariable>
  // {
  //   using type = OldStoredTemperatureVelocityVariable;
  // };

  template<typename Var> struct velocity_var {};

  // template<>
  // struct velocity_var<StoredVoltageVariable>
  // {
  //   using type = StoredElectricVelocityVariable;
  // };

  // template<>
  // struct velocity_var<OldStoredVoltageVariable>
  // {
  //   using type = OldStoredElectricVelocityVariable;
  // };
  
  // template<>
  // struct velocity_var<StoredTemperatureVariable>
  // {
  //   using type = StoredTemperatureVelocityVariable;
  // };

  // template<>
  // struct velocity_var<OldStoredTemperatureVariable>
  // {
  //   using type = OldStoredTemperatureVelocityVariable;
  // };
}

