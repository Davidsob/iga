#pragma once

#include "iga/constraints/ConstraintBase.h"
#include "iga/constraints/ConstraintTags.h"
#include "iga/operators/VectorizedVariable.h"
#include "iga/weakforms/WeakForms.h"

#include "shellmechanics/variables/StoredSixDofDisplacementVariable.h" 

#include "utils/MatrixTypes.h"

template<typename Var>
class DirichletPointConstraint 
  : public ConstraintBase
  , public PointConstraint 
{
  class Residual;
  class Jacobian;
public:

  ~DirichletPointConstraint() = default;

  WeakFormBase const * getResidual() const override
  { return dynamic_cast<WeakFormBase const *>(&_residual); }

  WeakFormBase const * getJacobian() const override
  { return dynamic_cast<WeakFormBase const *>(&_jacobian); }

  Residual const &getLinearForm() const { return _residual; }
  Jacobian const &getBilinearForm() const { return _jacobian; }

protected:

  explicit DirichletPointConstraint(double const x, size_t const dof, std::string const &name)
    : ConstraintBase(name)
    , _x(x), _dof(std::min(dof,Var::ndof))
    , _var()
    , _residual(*this)
    , _jacobian(*this)
  {}

  template<typename Index>
  size_t const gdof(Index const &i) const
  {
    return i.mapper.dof()[0]*Var::ndof + _dof;
  }

private:

  class Residual : public WeakForm<Residual, DynamicVectorR>
  {
  public:
    using parent_t = DirichletPointConstraint<Var>;
    using base_t = WeakForm<Residual, DynamicVectorR>;
    using value_t = StaticVectorR<3>;

    template<typename ...Args>
    Residual(parent_t const &parent)
      : base_t("DirichletPointConstraint_Residual")
      , _parent(parent)
    {}

    ~Residual() = default;

    template<typename Index>
    value_t const eval(Index const &i) const
    {
      auto const gdof = _parent.gdof(i);
      auto const var = _parent._var(i);
      auto const val = _parent._x - var[_parent._dof];
      return value_t({double(gdof),0,val});
    }

  private:
    parent_t const &_parent;
  };

  class Jacobian
    : public WeakForm<Jacobian, DynamicMatrixR>
  {
  public:
    using parent_t = DirichletPointConstraint<Var>;
    using base_t = WeakForm<Jacobian, DynamicMatrixR>;
    using value_t = StaticVectorR<3>;

    Jacobian(parent_t const &parent)
      : base_t("DirichletPointConstraint_Jacobian")
      , _parent(parent)
    {}

    ~Jacobian() = default;

    template<typename Index>
    value_t const eval(Index const &i) const
    {
      auto const gdof = _parent.gdof(i);
      return value_t(double(gdof),double(gdof),1);
    }

  private:
    parent_t const &_parent;
  };

  double const _x;
  size_t const _dof;
  VectorizedVariable<Var> const _var;
  Residual const _residual;
  Jacobian const _jacobian;
};

class PrescribedShellDisplacementU
  : public DirichletPointConstraint<SixDofDisplacementVariable>
{
  using base_t = DirichletPointConstraint<SixDofDisplacementVariable>;
public:

  explicit PrescribedShellDisplacementU(double const x)
    : base_t(x,0,"PrescribedShellDisplacementU")
  {}
  ~PrescribedShellDisplacementU() = default;
};

class PrescribedShellDisplacementV
  : public DirichletPointConstraint<SixDofDisplacementVariable>
{
  using base_t = DirichletPointConstraint<SixDofDisplacementVariable>;
public:
  
  explicit PrescribedShellDisplacementV(double const x)
    : base_t(x,1,"PrescribedShellDisplacementV")
  {}
  ~PrescribedShellDisplacementV() = default;
};

class PrescribedShellDisplacementW
  : public DirichletPointConstraint<SixDofDisplacementVariable>
{
  using base_t = DirichletPointConstraint<SixDofDisplacementVariable>;
public:
  
  explicit PrescribedShellDisplacementW(double const x)
    : base_t(x,2,"PrescribedShellDisplacementW")
  {}
  ~PrescribedShellDisplacementW() = default;
};
