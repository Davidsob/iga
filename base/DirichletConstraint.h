#pragma once

#include "ConstraintBase.h"
#include "MatrixTypes.h"
#include "WeakForms.h"

#include "StoredVoltageVariable.h"
#include "StoredTemperatureVariable.h"
#include "SimulationClock.h"

template<size_t N>
struct IdentityOperator
{
  using value_t = StaticMatrixR<N,N>;

  template<typename Index> 
  value_t const operator()(Index const &) const
  {
    return _eye;
  }

  static value_t const _eye;
};

template<size_t N>
StaticMatrixR<N,N> const IdentityOperator<N>::_eye = StaticMatrixR<N,N>::Identity();

template<typename T, typename Var>
class DirichletConstraint
  : public ConstraintBase
{
public:

private:
  class _f;

public:
  using linear_form = LinearWeakForm<IdentityOperator<1>, _f, DynamicVectorR>;
  using bilinear_form = BBWeakForm<IdentityOperator<1>, DynamicMatrixR>;

  template<typename ...Args>
  DirichletConstraint(Args ...args)
    : ConstraintBase("DirichletConstraint")
    , _load(args...)
    , _residual(IdentityOperator<1>(), _load)
    , _jacobian(IdentityOperator<1>())
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
    using value_t = StaticVectorR<1>;

    template<typename ...Args>
    _f(Args ...args) : _value(args...), _var() {}
    virtual ~_f() {};

    template<typename Index>
    value_t operator()(Index const &x) const
    {
      auto constraint = value_t(_value(x));
      auto variable = _var(x);
      auto diff = constraint - variable;
      // std::cout << "+++ Dirchlet op()" << std::endl;
      // std::cout << "constraint = " << constraint << std::endl;
      // std::cout << "variable = " << variable << std::endl;
      // std::cout << "delta = " << diff << std::endl;
      // std::cout << "--- Dirchlet op()\n" << std::endl;
      return diff;
      // return constraint;
    }

  private:
    T _value;
    Var _var;
  };

  _f const _load;
  linear_form const _residual;
  bilinear_form const _jacobian;
};


// template<typename T>
// class AppliedVoltage
//   : public DirichletConstraint<T, VoltageVariable>
// {
// public:
//   template<typename ...Args>
//   AppliedVoltage(Args ...args) : DirichletConstraint<T, VoltageVariable>(args...) {}
//   ~AppliedVoltage() = default;
// };

// template<typename T>
// class AppliedTemperature
//   : public DirichletConstraint<T, TemperatureVariable>
// {
// public:
//   template<typename ...Args>
//   AppliedTemperature(Args ...args) : DirichletConstraint<T, TemperatureVariable>(args...) {}
//   ~AppliedTemperature() = default;
// };