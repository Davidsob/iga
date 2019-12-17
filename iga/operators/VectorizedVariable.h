#pragma once

#include "utils/MatrixTypes.h"

template<typename Var>
class VectorizedVariable
{
public:

  using value_t = DynamicVectorR;

  VectorizedVariable()
    : _var() {}

  virtual ~VectorizedVariable() = default;

  template<typename Index>
  value_t const operator()(Index const &p) const
  {
    DynamicMatrixR const var = _var(p).transpose();
    Eigen::Map<DynamicVectorR const> const out(var.data(),var.size());
    return out;
  }

private:
  Var const _var;
};
