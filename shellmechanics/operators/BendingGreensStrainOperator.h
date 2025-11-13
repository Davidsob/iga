#pragma once

#include "iga/ManifoldElementMapper.h"
#include "iga/operators/Dot.h"
#include "iga/operators/VectorizedVariable.h"

#include "DisplacementIdOperator.h"
#include "RotationIdOperator.h"

#include "shellmechanics/variables/StoredSixDofDisplacementVariable.h"

#include "utils/MatrixTypes.h"

class BendingGreensStrainOperator 
{
public:
  using value_t = DynamicMatrixR;
  BendingGreensStrainOperator() : _var() {}; 
  ~BendingGreensStrainOperator() = default;

  template<typename Point>
  value_t const operator()(Point const &p) const
  {
    static DisplacementIdOperator const uid;
    static RotationIdOperator const rid;
    static Dot const dot;

    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&p.mapper);
    auto const grad = iga::CompactShapeFunctionDerivatives(p.para[0], p.para[1], mapper->manifold());
    auto const basis = mapper->covariantBasis(p.para);

    auto const dAd1 = rid(grad.row(0));
    auto const dAd2 = rid(grad.row(1));
    auto const dUd1 = uid(grad.row(0));
    auto const dUd2 = uid(grad.row(1));

    auto const var  = _var(p);
    auto const dud1 = dUd1*var;
    auto const dud2 = dUd2*var;
    auto const dad1 = dAd1*var;
    auto const dad2 = dAd2*var;

    StaticVectorR<3> const a1 = basis.col(0);
    StaticVectorR<3> const a2 = basis.col(1);

    auto const dof = grad.cols();
    DynamicMatrixR b(3,6*dof); b.setZero();
    b.row(0) = dot(a1,dAd1)+dot(dad1,dUd1)+dot(dud1,dAd1);
    b.row(1) = dot(a2,dAd2)+dot(dad2,dUd2)+dot(dud2,dAd2);
    b.row(2) = dot(a1,dAd2)+dot(a2,dAd1);
    b.row(2) += dot(dud1,dAd2) + dot(dud2,dAd1);
    b.row(2) += dot(dad1,dUd2) + dot(dad2,dUd1);

    return b;
  }

private:
  VectorizedVariable<SixDofDisplacementVariable> const _var;
};