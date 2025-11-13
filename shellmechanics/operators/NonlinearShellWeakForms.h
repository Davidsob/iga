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
#include "ShellStrain.h"
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
    auto const A = rid(shape);
    auto const n = mapper->covariantBasis(p.para).col(2);

    return dot(n,A);
  }
};

class NonlinearPenaltyStiffnessOperator
{
public:
  using value_t = DynamicRowVectorR;

  NonlinearPenaltyStiffnessOperator() : _var() {};
  virtual ~NonlinearPenaltyStiffnessOperator() = default;

  template<typename Index>
  value_t const operator()(Index const &p) const
  { 
    static RotationIdOperator const rid;
    static Dot const dot;
    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&p.mapper);
    auto const shape = p.mapper.shape(p.para);
    auto const A = rid(shape);

    auto const a = A*_var(p);
    auto const n = mapper->covariantBasis(p.para).col(2);
    auto const na = n+a;
    auto const mag = na.norm();

    return dot(na,A)/mag;
  }

  VectorizedVariable<SixDofDisplacementVariable> const _var;
};

template<typename MembraneTangent, typename ShearTangent>
class NonlinearMaterialStiffness
  : public WeakForm<NonlinearMaterialStiffness<MembraneTangent, ShearTangent>, DynamicMatrixR>
{
  using base_t          = WeakForm<NonlinearMaterialStiffness<MembraneTangent, ShearTangent>  , DynamicMatrixR>; 
  using membrane_form   = BdBWeakForm<MembraneGreensStrainOperator  , MembraneTangent      , DynamicMatrixR>;
  using bending_form    = BdBWeakForm<BendingGreensStrainOperator   , MembraneTangent      , DynamicMatrixR>;
  using transverse_form = BdBWeakForm<TransverseGreensStrainOperator, ShearTangent         , DynamicMatrixR>;
  using penalty_form    = BBWeakForm<LinearPenaltyStiffnessOperator, DynamicMatrixR>;
public:

  using value_t = DynamicMatrixR;

  NonlinearMaterialStiffness(MembraneTangent const &membrane,
                          ShearTangent const &shear,
                          double const thickness)
    : base_t("NonlinearMaterialStiffness")
    , _mf(MembraneGreensStrainOperator()  , membrane)
    , _bf(BendingGreensStrainOperator()   , membrane)
    , _tf(TransverseGreensStrainOperator(), shear)
    , _pf(LinearPenaltyStiffnessOperator())
    , _thickness(thickness)
  {}

  virtual ~NonlinearMaterialStiffness() = default;

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


template<typename MembraneTangent, typename ShearTangent>
class NonlinearMaterialResidual 
  : public WeakForm<NonlinearMaterialResidual<MembraneTangent, ShearTangent>, DynamicVectorR>
{
  using base_t          = WeakForm<NonlinearMaterialResidual<MembraneTangent, ShearTangent>  , DynamicVectorR>; 
  using membrane_form   = LinearWeakForm<MembraneGreensStrainOperator  , MembranePK2Stress       , DynamicVectorR>;
  using bending_form    = LinearWeakForm<BendingGreensStrainOperator   , BendingPK2Stress        , DynamicVectorR>;
  using transverse_form = LinearWeakForm<TransverseGreensStrainOperator, TransverseShearPK2Stress, DynamicVectorR>;

public:

  using value_t = DynamicVectorR;

  NonlinearMaterialResidual(MembraneTangent const &membrane,
                            ShearTangent const &shear,
                            double const thickness)
    : base_t("NonlinearMaterialResidual")
    , _msig(membrane, MembraneGreensStrain())
    , _bsig(membrane, BendingGreensStrain())
    , _tsig(shear,    TransverseShearGreensStrain())
    , _mf(MembraneGreensStrainOperator()  , _msig)
    , _bf(BendingGreensStrainOperator()   , _bsig)
    , _tf(TransverseGreensStrainOperator(), _tsig)
    , _thickness(thickness)
    , _bp()
    , _var()
  {}

  virtual ~NonlinearMaterialResidual() = default;

  template<typename Index>
  value_t const eval(Index const &i) const
  {
    double const penalty{1e8};
    static RotationIdOperator const rid;
    // compute the residual of the penalty term
    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&i.mapper);
    auto const shape = i.mapper.shape(i.para);
    auto const A = rid(shape);

    StaticVectorR<3> const a = A*_var(i);
    StaticVectorR<3> const n = mapper->covariantBasis(i.para).col(2);
    auto const f = n.dot(a);
    value_t const rp{f*_bp(i).transpose()};
    // compute the other residuals
    value_t const rm{_mf.eval(i)};
    value_t const rb{_bf.eval(i)};
    value_t const rt{_tf.eval(i)};
    // return the negative of the sum of the residals
    return -(_thickness*(rm+rt) + (std::pow(_thickness,3.0)/12.0)*rb + penalty*rp);
  }

private:
  MembranePK2Stress        const _msig;
  BendingPK2Stress         const _bsig;
  TransverseShearPK2Stress const _tsig;

  membrane_form   const _mf;
  bending_form    const _bf;
  transverse_form const _tf;
  double          const _thickness;

  LinearPenaltyStiffnessOperator              const _bp;
  VectorizedVariable<SixDofDisplacementVariable> const _var;
};

template<typename MembraneTangent, typename ShearTangent>
class NonlinearGeometricStiffness
  : public WeakForm<NonlinearGeometricStiffness<MembraneTangent, ShearTangent>, DynamicMatrixR>
{
  using MembraneStress  = MembranePK2Stress;
  using BendingStress    = BendingPK2Stress;
  using TransverseStress = TransverseShearPK2Stress;
  using base_t           = WeakForm<NonlinearGeometricStiffness<MembraneTangent, ShearTangent>  , DynamicMatrixR>; 
public:

  using value_t = DynamicMatrixR;

  NonlinearGeometricStiffness(MembraneTangent const &membrane,
                              ShearTangent const &shear,
                              double const thickness)
    : base_t("NonlinearGeometricStiffness")
    , _msig(membrane, MembraneGreensStrain())
    , _bsig(membrane, BendingGreensStrain())
    , _tsig(shear,    TransverseShearGreensStrain())
    , _thickness(thickness)
    // , _bp()
    // , _var()
  {}

  virtual ~NonlinearGeometricStiffness() = default;

  template<typename Index>
  value_t const eval(Index const &i) const
  {
    // double const penalty{1e8};
    value_t const km{_membrane(i)};
    value_t const kb{_bending(i)};
    value_t const kt{_transverse(i)};
    // value_t const kp{_penalty(i)};
    // return _thickness*(km+kt) + (std::pow(_thickness,3.0)/12.0)*kb + 0.0*penalty*kp;
    return _thickness*(km+kt) + (std::pow(0*_thickness,3.0)/12.0)*kb;
  }

private:

  template<typename Index>
  value_t const _membrane(Index const &i) const
  {
    static DisplacementIdOperator const uid;
    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&i.mapper);
    auto const grad = iga::CompactShapeFunctionDerivatives(i.para[0], i.para[1], mapper->manifold());
    auto const dUd1 = uid(grad.row(0));
    auto const dUd2 = uid(grad.row(1));
    auto const sig = _msig(i);

    return sig[0]*dUd1.transpose()*dUd1
         + sig[1]*dUd2.transpose()*dUd2
         + sig[2]*(dUd1.transpose()*dUd2 + dUd2.transpose()*dUd1);
  }

  template<typename Index>
  value_t const _bending(Index const &i) const
  {
    static DisplacementIdOperator const uid;
    static RotationIdOperator     const rid;
    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&i.mapper);
    auto const grad = iga::CompactShapeFunctionDerivatives(i.para[0], i.para[1], mapper->manifold());

    auto const dUd1 = uid(grad.row(0));
    auto const dUd2 = uid(grad.row(1));
    auto const dAd1 = rid(grad.row(0));
    auto const dAd2 = rid(grad.row(1));
    auto const sig = _bsig(i);

    return sig[0]*(dUd1.transpose()*dAd1 + dAd1.transpose()*dUd1)
         + sig[1]*(dUd2.transpose()*dAd2 + dAd2.transpose()*dUd2)
         + sig[2]*(dUd1.transpose()*dAd2 + dUd2.transpose()*dAd1)
         + sig[2]*(dAd1.transpose()*dUd2 + dAd2.transpose()*dUd1);
  }

  template<typename Index>
  value_t const _transverse(Index const &i) const
  {
    static DisplacementIdOperator const uid;
    static RotationIdOperator     const rid;
    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&i.mapper);
    auto const grad = iga::CompactShapeFunctionDerivatives(i.para[0], i.para[1], mapper->manifold());
    auto const shape = i.mapper.shape(i.para);

    auto const dUd1 = uid(grad.row(0));
    auto const dUd2 = uid(grad.row(1));
    auto const A    = rid(shape);
    auto const sig = _tsig(i);

    return sig[0]*(dUd2.transpose()*A + A.transpose()*dUd2)
         + sig[1]*(dUd1.transpose()*A + A.transpose()*dUd1);
  }

  // template<typename Index>
  // value_t const _penalty(Index const &i) const
  // {
  //   static RotationIdOperator     const rid;
  //   auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&i.mapper);
  //   auto const shape = i.mapper.shape(i.para);
  //   auto const A    = rid(shape);

  //   StaticVectorR<3> const a = A*_var(i);
  //   StaticVectorR<3> const n = mapper->covariantBasis(i.para).col(2);
  //   StaticVectorR<3> const an = a+n;
  //   auto const mag = an.norm();
  //   auto const f = mag-1.0;
  //   auto const g = _bp(i);

  //   return (f/mag)*(A.transpose()*A + g.transpose()*g);
  // }

  MembraneStress   const _msig;
  BendingStress    const _bsig;
  TransverseStress const _tsig;
  double           const _thickness;

  // NonlinearPenaltyStiffnessOperator const _bp;
  // VectorizedVariable<SixDofDisplacementVariable> const _var;
};
