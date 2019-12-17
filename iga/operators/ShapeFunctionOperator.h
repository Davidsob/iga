#pragma once

#include "utils/MatrixTypes.h"

class ShapeFunctionOperator 
{
public:
  using value_t = DynamicRowVectorR;

  ShapeFunctionOperator() = default;
  virtual ~ShapeFunctionOperator() = default;

  template<typename Index>
  value_t const operator()(Index const &p) const
  {
    return p.mapper.shape(p.para);
  }
};
