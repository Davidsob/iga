#pragma once

#include "utils/MatrixTypes.h"

class RotationIdOperator 
{
public:
  using value_t = DynamicMatrixR;

  RotationIdOperator() = default;
  virtual ~RotationIdOperator() = default;

  value_t const operator()(DynamicRowVectorR const &X) const
  {
    auto const dof = X.size();
    value_t op(value_t::Zero(3,6*dof));
    for (int i = 0; i < dof; i ++)
    {
      for (size_t j = 0; j < 3; j++) op(j,6*i+3+j) = X[i];
    }
    return op;
  }
};


