#pragma once

#include "utils/MatrixTypes.h"

class DisplacementIdOperator 
{
public:
  using value_t = DynamicMatrixR;
  
  DisplacementIdOperator() = default;
  virtual ~DisplacementIdOperator() = default;

  value_t const operator()(DynamicRowVectorR const &X) const
  {
    auto const dof = X.size();
    value_t op(value_t::Zero(3,6*dof));
    for (int i = 0; i < dof; i ++)
    {
      for (size_t j = 0; j < 3; j++) op(j,6*i+j) = X[i];
    }
    return op;
  }
};


