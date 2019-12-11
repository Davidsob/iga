#pragma once

#include "ConstraintBase.h"

#include "weakforms/WeakForms.h"

// #include "StoredVoltageVariable.h"
// #include "StoredTemperatureVariable.h"
#include "base/SimulationClock.h"

#include "utils/MatrixTypes.h"

#include <cassert>

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

template<size_t NDOF>
struct BIdOperator
{
  using value_t = DynamicMatrixR;

  template<typename Index>
  value_t const operator()(Index const &p) const
  {
    DynamicVectorR const shape = p.mapper.shape(p.para);
    value_t b(value_t::Zero(NDOF,NDOF*shape.rows()));
    for (int i = 0; i < shape.rows(); i++)
    {
      for (size_t j = 0; j < NDOF; j++)
      {
        b(j,NDOF*i+j) = shape[i];
      }
    }
    return b;
  }
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
  using id_op = BIdOperator<Var::ndof>;
  using linear_form = LinearWeakForm<id_op, _f, DynamicVectorR>;
  using bilinear_form = BBWeakForm<id_op, DynamicMatrixR>;

  template<typename ...Args>
  DirichletConstraint(Args ...args)
    : ConstraintBase("DirichletConstraint")
    , _load(args...)
    , _residual(id_op(), _load)
    , _jacobian(id_op())
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
    value_t operator()(Index const &p) const
    {
      // auto const shape = p.mapper.shapeFunctions(p.para);
      auto constraint = value_t(_value(p));
      auto variable = value_t(0);
      // auto variable = _var(x);
      value_t diff = constraint - variable;
      // std::cout << "+++ Dirchlet op()" << std::endl;
      // std::cout << "constraint = " << constraint << std::endl;
      // std::cout << "variable = " << variable << std::endl;
      // std::cout << "delta = " << diff << std::endl;
      // std::cout << "--- Dirchlet op()\n" << std::endl;
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