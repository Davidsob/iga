#pragma once

#include "utils/MatrixTypes.h"

template<size_t NDOF>
struct BIdOperator
{
  using value_t = SparseMatrixR;

  template<typename Index>
  value_t const operator()(Index const &p) const
  {
    DynamicVectorR const shape = p.mapper.shape(p.para);
    DynamicMatrixR b(DynamicMatrixR::Zero(NDOF,NDOF*shape.rows()));
    for (int i = 0; i < shape.rows(); i++)
    {
      for (size_t j = 0; j < NDOF; j++)
      {
        b(j,NDOF*i+j) = shape[i];
      }
    }
    return b.sparseView();
  }
};




