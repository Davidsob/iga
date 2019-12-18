#pragma once

#include "iga/ManifoldElementMapper.h"
#include "iga/operators/Dot.h"

#include "DisplacementIdOperator.h"
#include "RotationIdOperator.h"

#include "utils/MatrixTypes.h"

class LinearTransverseStrainOperator 
{
public:
  using value_t = DynamicMatrixR;
  LinearTransverseStrainOperator() = default; 
  ~LinearTransverseStrainOperator() = default;

  template<typename Point>
  value_t const operator()(Point const &p) const
  {
    static DisplacementIdOperator const uid;
    static RotationIdOperator const rid;
    static Dot const dot;

    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&p.mapper);
    auto const shape = iga::CompactShapeFunctions(p.para[0], p.para[1], mapper->manifold());
    auto const grad  = iga::CompactShapeFunctionDerivatives(p.para[0], p.para[1], mapper->manifold());
    auto const basis = mapper->covariantBasis(p.para);

    auto const dNd1 = uid(grad.row(0));
    auto const dNd2 = uid(grad.row(1));
    auto const N    = rid(shape);

    StaticVectorR<3> const a1 = basis.col(0);
    StaticVectorR<3> const a2 = basis.col(1);
    StaticVectorR<3> const n  = basis.col(2);

    auto const dof = grad.cols();
    DynamicMatrixR b(2,6*dof); b.setZero();

    b.row(0) = dot(n,dNd2) + dot(a2,N);
    b.row(1) = dot(n,dNd1) + dot(a1,N);

    return b;
  }
};