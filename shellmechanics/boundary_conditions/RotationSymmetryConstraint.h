#pragma once

#include "iga/CurveElementMapper.h"
#include "iga/constraints/ConstraintBase.h"
#include "iga/constraints/ConstraintTags.h"
#include "iga/operators/BIdOperator.h"
#include "iga/operators/Dot.h"
#include "iga/weakforms/WeakForms.h"

#include "utils/MatrixTypes.h"

#include "shellmechanics/operators/DisplacementIdOperator.h"
#include "shellmechanics/operators/RotationIdOperator.h"

#include "NormalSourceBase.h"

class RotationSymmetryConstraint 
  : public ConstraintBase
  , public ScalarConstraint
{
  class Residual;
  class Jacobian;
public:

  template<typename ...Args>
  RotationSymmetryConstraint(Args const &...args)
    : ConstraintBase("ShellRotationSymmetryConstraint")
    , _residual(args...)
    , _jacobian(args...)
  {}

  ~RotationSymmetryConstraint() = default;

  WeakFormBase const * getResidual() const override
  { return dynamic_cast<WeakFormBase const *>(&_residual); }

  WeakFormBase const * getJacobian() const override
  { return dynamic_cast<WeakFormBase const *>(&_jacobian); }

  Residual const &getLinearForm() const { return _residual; }
  Jacobian const &getBilinearForm() const { return _jacobian; }


private:

  class Residual : public WeakForm<Residual, DynamicVectorR>
  {
  public:
    using value_t = DynamicVectorR;
    using base_t = WeakForm<Residual, DynamicVectorR>;

    template<typename ...Args>
    Residual(Args &&...args)
      : base_t("ShellRotationSymmetryConstraint_Residual")
    {}

    ~Residual() = default;

    template<typename Index>
    value_t const eval(Index const &i) const
    {
      return DynamicVectorR::Zero(1);
    }

  };

  class Jacobian
    : public WeakForm<Jacobian, DynamicMatrixR>
  {
  public:
    using value_t = DynamicMatrixR;
    using base_t = WeakForm<Jacobian, DynamicMatrixR>;

    Jacobian(double uv, size_t dir, NormalSourceBase const &normal)
      : base_t("ShellRotationSymmetryConstraint_Jacobian")
      , _uv(uv), _direction(dir), _normal(normal)
    {}

    ~Jacobian() = default;

    template<typename Index>
    value_t const eval(Index const &i) const
    {
      static RotationIdOperator const rid;
      static Dot const dot;

      auto              const mapper = dynamic_cast<CurveElementMapper const *>(&i.mapper);
      DynamicRowVectorR const shape  = mapper->shape(i.para);
      StaticVectorR<3>  const t1 = mapper->tangent(i.para);
      StaticVectorR<3>  const n  = _direction ? _normal(_uv,i.para[0]) : _normal(i.para[0],_uv);
      StaticVectorR<3>  const t2 = n.cross(t1);

      value_t const na{dot(t2,rid(shape))};

      return na;
    }

  private:
    double _uv;
    size_t _direction;
    NormalSourceBase const &_normal;
  };

  Residual const _residual;
  Jacobian const _jacobian;
};
