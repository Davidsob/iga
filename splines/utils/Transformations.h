#pragma once

#include "Quaternion.h"
#include "VectorOperations.h"
#include <iostream>

using namespace vector_ops;

namespace transform
{
  template<typename T>
  T translate(T const &original, std::vector<double> const &translation)
  {
    T copy(original);
    for (auto &x : copy.Q) x += translation;
    return copy;
  }

  template<typename T>
  T rotate(T const &original,
          std::vector<double> const &n0,
          std::vector<double> const &origin,
          std::vector<double> const &n)
  {
    // make quaternion for roatation
    auto C = dot(n0,n)/(norm(n0)*norm(n));
    auto th = -std::acos(C);
    auto qv = cross(n0,n); qv = normalize(qv);
    Quaternion q(std::cos(0.5*th), std::sin(0.5*th)*qv);

    T copy(original);

    for (auto &x : copy.Q)
    {
      x = rotate(q,Quaternion(0,x)).vectorComponent() + origin;
    }
    return copy;
  }

  template<typename T>
  T project(T const &original,
          std::vector<double> const &N0,
          std::vector<double> const &origin,
          std::vector<double> &n)
  {
    // make quaternion for roatation
    T copy(original);
    n = normalize(n);
    for (auto &x : copy.Q)
    {
      auto u = origin-x;
      auto s = dot(n,u)/dot(n,N0);
      x += s*N0;
    }
    return copy;
  }
}