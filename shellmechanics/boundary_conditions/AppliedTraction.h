#pragma once

#include "base/Values.h"

#include "iga/operators/ShapeFunctionOperator.h"

#include "shellmechanics/operators/DisplacementIdOperator.h"

#include "utils/MatrixTypes.h"

template<typename Traction>
class AppliedTraction 
  : public BoundaryConditionBase 
{
  using BId_t = MappedValue<ShapeFunctionOperator, DisplacementIdOperator, DynamicMatrixR>;
  using linear_form = LinearWeakForm<BId_t, Traction, DynamicVectorR>;

public:
  template<typename ...Args>
  AppliedTraction(Args ...args)
    : BoundaryConditionBase("AppliedTraction")
    , _traction(args...)
    , _traction_form(BId_t(), _traction)
  {}

  virtual ~AppliedTraction() = default;

  WeakFormBase const * getRhs() const override
  { return dynamic_cast<WeakFormBase const *>(&_traction_form); }

private:

  Traction    const _traction;
  linear_form const _traction_form;
};
