#pragma once

#include "iga/constraints/DirichletConstraint.h"

#include "utils/MatrixTypes.h"


class FixedFixed
  : public DirichletConstraint<ConstantValue<StaticVectorR<6>>, SixDofDisplacementVariable>
{
  using base_t = DirichletConstraint<ConstantValue<StaticVectorR<6>>, SixDofDisplacementVariable>;
public:
  using value_t = StaticVectorR<6>;
  FixedFixed()
    : base_t(value_t::Zero())
  {}
  ~FixedFixed() = default;
};
