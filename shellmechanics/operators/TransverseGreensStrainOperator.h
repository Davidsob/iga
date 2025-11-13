#pragma once

#include "iga/ManifoldElementMapper.h"
#include "iga/operators/Dot.h"
#include "iga/operators/VectorizedVariable.h"

#include "DisplacementIdOperator.h"
#include "RotationIdOperator.h"

#include "shellmechanics/variables/StoredSixDofDisplacementVariable.h"

#include "utils/MatrixTypes.h"

class TransverseGreensStrainOperator 
{
public:
  using value_t = DynamicMatrixR;
  TransverseGreensStrainOperator() = default; 
  ~TransverseGreensStrainOperator() = default;

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

    auto const dUd1 = uid(grad.row(0));
    auto const dUd2 = uid(grad.row(1));
    auto const A    = rid(shape);

    auto const var  = _var(p);
    auto const dud1 = dUd1*var;
    auto const dud2 = dUd2*var;
    auto const a    = A*var;

    StaticVectorR<3> const a1 = basis.col(0);
    StaticVectorR<3> const a2 = basis.col(1);
    StaticVectorR<3> const n  = basis.col(2);

    auto const dof = grad.cols();
    DynamicMatrixR b(2,6*dof); b.setZero();

    b.row(0) = dot(n,dUd2) + dot(a2,A) + dot(a,dUd2) + dot(dud2,A);
    b.row(1) = dot(n,dUd1) + dot(a1,A) + dot(a,dUd1) + dot(dud1,A);

    return b;
  }

private:
  VectorizedVariable<SixDofDisplacementVariable> const _var;
};