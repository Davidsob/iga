#pragma once

#include "iga/ManifoldElementMapper.h"
#include "iga/operators/Dot.h"

#include "RotationIdOperator.h"

#include "utils/MatrixTypes.h"

class LinearBendingStrainOperator 
{
public:
  using value_t = DynamicMatrixR;
  LinearBendingStrainOperator() = default; 
  ~LinearBendingStrainOperator() = default;

  template<typename Point>
  value_t operator()(Point const &p) const
  {
    static RotationIdOperator const rid;
    static Dot const dot;

    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&p.mapper);
    auto const grad = iga::CompactShapeFunctionDerivatives(p.para[0], p.para[1], mapper->manifold());
    auto const basis = mapper->covariantBasis(p.para);

    auto const dNd1 = rid(grad.row(0));
    auto const dNd2 = rid(grad.row(1));

    StaticVectorR<3> const a1 = basis.col(0);
    StaticVectorR<3> const a2 = basis.col(1);

    auto const dof = grad.cols();
    DynamicMatrixR b(3,6*dof); b.setZero();
    b.row(0) = dot(a1,dNd1);
    b.row(1) = dot(a2,dNd2);
    b.row(2) = dot(a1,dNd2) + dot(a2,dNd1);

    return b;
  }
};