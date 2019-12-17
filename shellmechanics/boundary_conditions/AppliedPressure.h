#pragma once

#include "base/Values.h"

#include "iga/operators/ShapeFunctionOperator.h"

#include "shellmechanics/operators/DisplacementIdOperator.h"

#include "utils/MatrixTypes.h"

template<typename Pressure>
class AppliedPressure 
  : public BoundaryConditionBase 
{
  class PressureOperator;
  using BId_t = FunctorMappedValue<ShapeFunctionOperator, DisplacementIdOperator, DynamicMatrixR>;
  using linear_form = LinearWeakForm<BId_t, PressureOperator, DynamicVectorR>;

public:
  template<typename ...Args>
  AppliedPressure(Args ...args)
    : BoundaryConditionBase("AppliedPressure")
    , _pressure(args...)
    , _pressure_form(BId_t(), _pressure)
  {}

  virtual ~AppliedPressure() = default;

  WeakFormBase const * getRhs() const override
  { return dynamic_cast<WeakFormBase const *>(&_pressure_form); }

private:
  class PressureOperator
  {
  public:
    using value_t = StaticVectorR<3>;

    template<typename ...Args>
    PressureOperator(Args ...args) : _pressure(args...) {}
    virtual ~PressureOperator() {};

    template<typename Index>
    value_t operator()(Index const &i) const
    {
      auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&i.mapper);
      value_t const n = mapper->covariantBasis(i.para).col(2);
      return -_pressure(i)*n;
    }

  private:
    Pressure _pressure;
  };

  PressureOperator const _pressure;         
  linear_form const _pressure_form;
};
