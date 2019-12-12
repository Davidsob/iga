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

  RampedValue(T const x, double s=1.0, double b=0.0)
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
  double const _b;
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
