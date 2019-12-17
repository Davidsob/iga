#pragma once

#include "base/Values.h"
#include "iga/operators/VectorizedVariable.h"
#include "iga/weakforms/WeakForms.h"

#include "utils/MatrixTypes.h"

#include "iga/operators/Dot.h"

#include "RotationIdOperator.h"
#include "LinearMembraneStrainOperator.h"
#include "LinearBendingStrainOperator.h"
#include "LinearTransverseStrainOperator.h"


class LinearPenaltyStiffnessOperator
{
public:
  using value_t = DynamicRowVectorR;

  LinearPenaltyStiffnessOperator() {};
  virtual ~LinearPenaltyStiffnessOperator() = default;

  template<typename Index>
  value_t const operator()(Index const &p) const
  { 
    static RotationIdOperator const rid;
    static Dot const dot;
    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&p.mapper);
    auto const shape = p.mapper.shape(p.para);
    auto const n = mapper->covariantBasis(p.para).col(2);
    return dot(n,rid(shape));
  }
};


template<typename MembraneTangent, typename ShearTangent>
class LinearMaterialStiffness
  : public WeakForm<LinearMaterialStiffness<MembraneTangent, ShearTangent>, DynamicMatrixR>
{
  using base_t          = WeakForm<LinearMaterialStiffness<MembraneTangent, ShearTangent>  , DynamicMatrixR>; 
  using membrane_form   = BdBWeakForm<LinearMembraneStrainOperator  , MembraneTangent      , DynamicMatrixR>;
  using bending_form    = BdBWeakForm<LinearBendingStrainOperator   , MembraneTangent      , DynamicMatrixR>;
  using transverse_form = BdBWeakForm<LinearTransverseStrainOperator, ShearTangent         , DynamicMatrixR>;
  using penalty_form    = BBWeakForm<LinearPenaltyStiffnessOperator, DynamicMatrixR>;
public:

  using value_t = DynamicMatrixR;

  LinearMaterialStiffness(MembraneTangent const &membrane,
                          ShearTangent const &shear,
                          double const thickness)
    : base_t("LinearMaterialStiffness")
    , _mf(LinearMembraneStrainOperator()  , membrane)
    , _bf(LinearBendingStrainOperator()   , membrane)
    , _tf(LinearTransverseStrainOperator(), shear)
    , _pf(LinearPenaltyStiffnessOperator())
    , _thickness(thickness)
  {}

  virtual ~LinearMaterialStiffness() = default;

  template<typename Index>
  value_t const eval(Index const &i) const
  {
    double const penalty{1e8};
    value_t const km{_mf.eval(i)};
    value_t const kb{_bf.eval(i)};
    value_t const kt{_tf.eval(i)};
    value_t const kp{_pf.eval(i)};
    return _thickness*(km+kt) + (std::pow(_thickness,3.0)/12.0)*kb + penalty*kp;
  }

private:
  membrane_form   const _mf;
  bending_form    const _bf;
  transverse_form const _tf;
  penalty_form    const _pf;
  double          const _thickness;
};