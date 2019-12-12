#pragma once

#include "thermal/StoredTemperatureVariable.h"
#include "thermal/StoredTemperatureRateVariable.h"

namespace StoredVariableTraits
{
  // class StoredTemperatureVariable;
  // class StoredTemperatureVariableIncrement;
  // class OldStoredTemperatureVariable;

  // class StoredTemperatureRateVariable;
  // class StoredTemperatureRateVariableIncrement;
  // class OldStoredTemperatureRateVariable;

  template<typename Var> struct increment_var {};

  template<>
  struct increment_var<StoredTemperatureVariable>
  {
    using type = StoredTemperatureVariableIncrement;
  };

  template<>
  struct increment_var<StoredTemperatureRateVariable>
  {
    using type = StoredTemperatureRateVariableIncrement;
  };

  template<typename Var> struct old_var {};

  template<>
  struct old_var<StoredTemperatureVariable>
  {
    using type = OldStoredTemperatureVariable;
  };

  template<>
  struct old_var<StoredTemperatureRateVariable>
  {
    using type = OldStoredTemperatureRateVariable;
  };

  template<typename Var> struct velocity_var {};

  template<>
  struct velocity_var<StoredTemperatureVariable>
  {
    using type = StoredTemperatureRateVariable;
  };

  template<>
  struct velocity_var<OldStoredTemperatureVariable>
  {
    using type = OldStoredTemperatureRateVariable;
  };
}

