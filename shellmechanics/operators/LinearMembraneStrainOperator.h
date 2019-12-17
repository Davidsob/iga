#pragma once

#include "iga/ManifoldElementMapper.h"
#include "iga/operators/Dot.h"

#include "DisplacementIdOperator.h"

#include "utils/MatrixTypes.h"

class LinearMembraneStrainOperator 
{
public:
  using value_t = DynamicMatrixR;

  LinearMembraneStrainOperator() = default; 
  ~LinearMembraneStrainOperator() = default;

  template<typename Point>
  value_t operator()(Point const &p) const
  {
    static DisplacementIdOperator const uid;
    static Dot const dot;

    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&p.mapper);
    auto const grad = iga::CompactShapeFunctionDerivatives(p.para[0], p.para[1], mapper->manifold());
    auto const basis = mapper->covariantBasis(p.para);

    auto const dNd1 = uid(grad.row(0));
    auto const dNd2 = uid(grad.row(1));

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