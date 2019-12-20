#pragma once

#include "base/Values.h"

#include "iga/operators/IdentityOperator.h"
#include "iga/operators/ShapeFunctionOperator.h"

#include "shellmechanics/operators/DisplacementIdOperator.h"

#include "utils/MatrixTypes.h"

template<typename Load>
class PointLoad 
  : public BoundaryConditionBase 
{
  using BId_t = MappedValue<IdentityOperator<1>, DisplacementIdOperator, DynamicMatrixR>;
  using linear_form = LinearWeakForm<BId_t, Load, DynamicVectorR>;

public:
  template<typename ...Args>
  PointLoad(Args ...args)
    : BoundaryConditionBase("PointLoad")
    , _load(args...)
    , _load_form(BId_t(), _load)
  {}

  virtual ~PointLoad() = default;

  WeakFormBase const * getRhs() const override
  { return dynamic_cast<WeakFormBase const *>(&_load_form); }

private:

  Load        const _load;
  linear_form const _load_form;
};
