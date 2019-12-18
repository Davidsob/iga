#pragma once

#include "utils/MatrixTypes.h"

#include <cmath>

template<typename T>
class ConstantValue
{
public:
  using value_t = T;

  ConstantValue(T const x) : _x(x) {}
  virtual ~ConstantValue() = default;

  template<typename Index>
  T operator()(Index const &i) const { return _x; }

private:
  T const _x;
};

template<typename T>
class NegativeValue
{
public:
  using value_t = typename T::value_t;

  template<typename ...Args>
  NegativeValue(Args ...args) : _x(args...) {}

  virtual ~NegativeValue() = default;

  template<typename Index>
  value_t operator()(Index const &i) const
  { 
    return -_x(i);
  }

private:
  T const _x;
};

template<typename A, typename B, typename Operator, typename Return>
class BinaryMappedValue 
{
public:
  using value_t = Return;

  BinaryMappedValue()
    : _a(), _b(), _op() {}

  BinaryMappedValue(A const a, B const b)
    : _a(a), _b(b), _op() {}

  virtual ~BinaryMappedValue() = default;

  template<typename Index>
  value_t const operator()(Index const &i) const
  { 
    return _op(_a(i),_b(i));
  }

private:
  A const _a;
  B const _b;
  Operator const _op;
};

template<typename A, typename Operator, typename Return>
class MappedValue 
{
public:
  using value_t = Return;

  MappedValue()
    : _a(), _op(){}

  MappedValue(A const &a)
    : _a(a), _op() {}

  virtual ~MappedValue() = default;

  template<typename Index>
  value_t operator()(Index const &i) const
  { 
    return _op(_a(i));
  }

private:
  A const _a;
  Operator const _op;
};

template<typename T>
class ScalarScaledValue 
{
public:
  using value_t = typename T::value_t;

  template<typename ...Args>
  ScalarScaledValue(double scale, Args ...args)
    : _scale(scale), _x(args...) {}

  virtual ~ScalarScaledValue() = default;

  template<typename Index>
  value_t operator()(Index const &i) const
  { 
    return _scale*_x(i);
  }

private:
  double const _scale;
  T const _x;
};

template<typename T>
class TransposedValue
{
public:
  using value_t = DynamicMatrixR;

  template<typename ...Args>
  TransposedValue(Args ...args) : _x(args...) {}

  virtual ~TransposedValue() = default;

  template<typename Index>
  value_t operator()(Index const &i) const
  { 
    return _x(i).transpose();
  }

private:
  T const _x;
};

template<typename T>
class RampedValue
{
public:
  using value_t = T;

  RampedValue(T const x, double s=1.0, T b=0.0)
    : _x(x), _s(s), _b(b)
  {}

  virtual ~RampedValue() = default;

  template<typename Index>
  T operator()(Index const &i) const
  {
    return (_s*i.time)*_x + _b; 
  }

private:
  T const _x;
  double const _s;
  T const _b;
};

template<typename T>
class SteppedValue
{
public:
  using value_t = T;

  SteppedValue(T const x, double rise_time, double slope)
    : _x(x), _rise_time(rise_time), _slope(slope)
  {}

  virtual ~SteppedValue() = default;

  template<typename Index>
  T operator()(Index const &i) const
  {
    return 0.5*(std::tanh(_slope*(i.time - _rise_time)) + 1.0)*_x; 
  }

private:
  T const _x;
  double const _rise_time;
  double const _slope;
};

using ConstantScalarValue = ConstantValue<double>;
using RampedScalarValue = RampedValue<double>;
using SteppedScalarValue = SteppedValue<double>;
