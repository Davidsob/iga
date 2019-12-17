#pragma once

#include "ConstraintBase.h"

#include "weakforms/WeakForms.h"

// #include "StoredVoltageVariable.h"
// #include "StoredTemperatureVariable.h"
#include "base/SimulationClock.h"

#include "utils/MatrixTypes.h"

#include <cassert>

#include "operators/BIdOperator.h"

template<typename T, typename Var, typename Bid = BIdOperator<Var::ndof>>
class DirichletConstraint
  : public ConstraintBase
{
public:

private:
  class _f;

public:
  using linear_form = LinearWeakForm<Bid, _f, DynamicVectorR>;
  using bilinear_form = BBWeakForm<Bid, DynamicMatrixR>;

  template<typename ...Args>
  DirichletConstraint(Args ...args)
    : ConstraintBase("DirichletConstraint")
    , _load(args...)
    , _residual(Bid(), _load)
    , _jacobian(Bid())
  {}

  virtual ~DirichletConstraint() = default;

  WeakFormBase const * getResidual() const override
  { return dynamic_cast<WeakFormBase const *>(&_residual); }

  WeakFormBase const * getJacobian() const override
  { return dynamic_cast<WeakFormBase const *>(&_jacobian); }

  linear_form const &getLinearForm() const { return _residual; }
  bilinear_form const &getBilinearForm() const { return _jacobian; }

private:

  class _f
  {
  public:
    using value_t = typename T::value_t;

    template<typename ...Args>
    _f(Args ...args) : _value(args...), _var() {}
    virtual ~_f() {};

    template<typename Index>
    value_t interpolatedVar(Index const &p) const
    {
      auto const shape = p.mapper.shape(p.para);
      return shape*_var(p);
    }

    template<typename Index>
    value_t operator()(Index const &p) const
    {
      auto constraint = value_t(_value(p));
      auto variable = interpolatedVar(p);
      auto diff = constraint - variable;
      return diff;
    }

  private:
    T _value;
    Var _var;
  };

  _f const _load;
  linear_form const _residual;
  bilinear_form const _jacobian;
};