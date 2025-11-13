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
#include "MembraneGreensStrainOperator.h"
#include "BendingGreensStrainOperator.h"
#include "TransverseGreensStrainOperator.h"
#include "ShellStress.h"

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
    // return _thickness*(km+kt) + (std::pow(_thickness,3.0)/12.0)*kb;
  }

private:
  membrane_form   const _mf;
  bending_form    const _bf;
  transverse_form const _tf;
  penalty_form    const _pf;
  double          const _thickness;
};

template<typename MembraneTangent, typename ShearTangent>
class MaterialResidual 
  : public WeakForm<MaterialResidual<MembraneTangent, ShearTangent>, DynamicVectorR>
{
  using base_t          = WeakForm<MaterialResidual<MembraneTangent, ShearTangent>  , DynamicVectorR>; 
  using membrane_form   = LinearWeakForm<LinearMembraneStrainOperator  , MembraneStress       , DynamicVectorR>;
  using bending_form    = LinearWeakForm<LinearBendingStrainOperator   , BendingStress        , DynamicVectorR>;
  using transverse_form = LinearWeakForm<LinearTransverseStrainOperator, TransverseShearStress, DynamicVectorR>;

public:

  using value_t = DynamicVectorR;

  MaterialResidual(MembraneTangent const &membrane,
                   ShearTangent const &shear,
                   double const thickness)
    : base_t("MaterialResidual")
    , _msig(MembraneStress(membrane, LinearMembraneStrain()))
    , _bsig(BendingStress(membrane , LinearBendingStrain()))
    , _tsig(TransverseShearStress(shear, LinearTransverseShearStrain()))
    , _mf(LinearMembraneStrainOperator()  , _msig)
    , _bf(LinearBendingStrainOperator()   , _bsig)
    , _tf(LinearTransverseStrainOperator(), _tsig)
    , _bp()
    , _var()
    , _thickness(thickness)
  {}

  virtual ~MaterialResidual() = default;

  template<typename Index>
  value_t const eval(Index const &i) const
  {
    // double const penalty{1e8};
    static RotationIdOperator const rid;
    // compute the other residuals
    value_t const rm{_mf.eval(i)};
    value_t const rb{_bf.eval(i)};
    value_t const rt{_tf.eval(i)};
    // compute the residual of the penalty term
    // auto mapper = dynamic_cast<ManifoldElementMapper const *>(&i.mapper);
    // StaticVectorR<3> const n = mapper->covariantBasis(i.para).col(2);
    // StaticVectorR<3> const a = rid(mapper->shape(i.para))*_var(i);
    // value_t const rp{n.dot(a)*_bp(i).transpose()};
    // return the negative of the sum of the residals
    // return -(_thickness*(rm+rt) + (std::pow(_thickness,3.0)/12.0)*rb + penalty*rp);
    return -(_thickness*(rm+rt) + (std::pow(_thickness,3.0)/12.0)*rb);
  }

private:
  MembraneStress          const _msig;
  BendingStress           const _bsig;
  TransverseShearStress   const _tsig;

  membrane_form   const _mf;
  bending_form    const _bf;
  transverse_form const _tf;
  LinearPenaltyStiffnessOperator _bp;
  VectorizedVariable<SixDofDisplacementVariable> const _var;
  double          const _thickness;
};
