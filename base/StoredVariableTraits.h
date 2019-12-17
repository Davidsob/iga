#pragma once

#include "thermal/StoredTemperatureVariable.h"
#include "thermal/StoredTemperatureRateVariable.h"
#include "solidmechanics/StoredDisplacementVariable.h"
#include "solidmechanics/StoredVelocityVariable.h"
#include "shellmechanics/variables/StoredSixDofDisplacementVariable.h"
#include "shellmechanics/variables/StoredSixDofVelocityVariable.h"

namespace StoredVariableTraits
{
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

  template<>
  struct increment_var<StoredDisplacementVariable<3>>
  {
    using type = StoredDisplacementVariableIncrement<3>;
  };

  template<>
  struct increment_var<StoredVelocityVariable<3>>
  {
    using type = StoredVelocityVariableIncrement<3>;
  };

  template<>
  struct increment_var<StoredSixDofDisplacementVariable>
  {
    using type = StoredSixDofDisplacementVariableIncrement;
  };

  template<>
  struct increment_var<StoredSixDofVelocityVariable>
  {
    using type = StoredSixDofVelocityVariableIncrement;
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

  template<>
  struct old_var<StoredDisplacementVariable<3>>
  {
    using type = OldStoredDisplacementVariable<3>;
  };

  template<>
  struct old_var<StoredVelocityVariable<3>>
  {
    using type = OldStoredVelocityVariable<3>;
  };

  template<>
  struct old_var<StoredSixDofDisplacementVariable>
  {
    using type = OldStoredSixDofDisplacementVariable;
  };

  template<>
  struct old_var<StoredSixDofVelocityVariable>
  {
    using type = OldStoredSixDofVelocityVariable;
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

  template<>
  struct velocity_var<StoredDisplacementVariable<3>>
  {
    using type = StoredVelocityVariable<3>;
  };

  template<>
  struct velocity_var<OldStoredDisplacementVariable<3>>
  {
    using type = OldStoredVelocityVariable<3>;
  };

  template<>
  struct velocity_var<StoredSixDofDisplacementVariable>
  {
    using type = StoredSixDofVelocityVariable;
  };

  template<>
  struct velocity_var<OldStoredSixDofDisplacementVariable>
  {
    using type = OldStoredSixDofVelocityVariable;
  };
}

