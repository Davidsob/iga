#pragma once

#include "utils/MatrixTypes.h"

template<size_t N>
struct IdentityOperator
{
  using value_t = StaticMatrixR<N,N>;

  template<typename Index> 
  value_t const operator()(Index const &) const
  {
    return _eye;
  }

  static value_t const _eye;
};

template<size_t N>
StaticMatrixR<N,N> const IdentityOperator<N>::_eye = StaticMatrixR<N,N>::Identity();


