#pragma once

#include <Eigen/Dense>
#include <type_traits>

#include <vector>

#include "utils/MatrixTypes.h"

namespace quadrature
{

  template<typename IntegrationPoint, typename Operator>
  inline void gauss(std::vector<IntegrationPoint> const &ips, std::vector<Triplet> &out, Operator && op)
  {
    for (auto const &ip : ips)
    {
      StaticVectorR<3> const x(ip.weight*ip.jdet*op(ip));
      out.push_back(convert::to<Triplet>(x));
    }
  }

  template<typename IntegrationPoint, typename Output, typename Operator>
  inline void gauss(std::vector<IntegrationPoint> const &ips, Output &out, Operator && op)
  {
    for (auto const &ip : ips)
    {
      out += ip.weight*ip.jdet*op(ip);
    }
  }
}