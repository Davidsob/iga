#pragma once

#include "iga/ManifoldElementMapper.h"
#include "iga/constraints/ConstraintBase.h"
#include "iga/constraints/ConstraintTags.h"
#include "iga/operators/Dot.h"
#include "iga/operators/VectorizedVariable.h"
#include "iga/weakforms/WeakForms.h"

#include "utils/MatrixTypes.h"

#include "shellmechanics/operators/RotationIdOperator.h"
#include "shellmechanics/variables/StoredSixDofDisplacementVariable.h"

class LinearShellConstraint 
  : public ConstraintBase
  , public ScalarConstraint
{
  class Residual;
  class Jacobian;
public:

  LinearShellConstraint()
    : ConstraintBase("LinearShellConstraint")
    , _residual()
    , _jacobian()
  {}

  ~LinearShellConstraint() = default;

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

    Residual()
      : base_t("LinearShellConstraint_Residual")
      , _var()
    {}

    ~Residual() = default;

    template<typename Index>
    value_t const eval(Index const &i) const
    {
      static RotationIdOperator const rid;
      auto mapper = dynamic_cast<ManifoldElementMapper const *>(&i.mapper);
      auto const shape = mapper->shape(i.para);
      auto const n = mapper->normal(i.para);
      auto const a = rid(shape)*_var(i);
      DynamicVectorR na(1); na[0] = n.dot(a);
      return -na;
    }

    VectorizedVariable<SixDofDisplacementVariable> const _var;

  };

  class Jacobian
    : public WeakForm<Jacobian, DynamicMatrixR>
  {
  public:
    using value_t = DynamicMatrixR;
    using base_t = WeakForm<Jacobian, DynamicMatrixR>;

    Jacobian()
      : base_t("LinearShellConstraint_Jacobian")
    {}

    ~Jacobian() = default;

    template<typename Index>
    value_t const eval(Index const &i) const
    {
      static RotationIdOperator const rid;
      static Dot const dot;

      auto mapper = dynamic_cast<ManifoldElementMapper const *>(&i.mapper);
      auto const shape = mapper->shape(i.para);
      auto const n = mapper->normal(i.para);
      auto const out = dot(n,rid(shape));
      return out;
    }
  };

  Residual const _residual;
  Jacobian const _jacobian;
};
