#pragma once

#include <Eigen/Dense>

#include <type_traits>

namespace quadrature
{

  template<typename IntegrationPoint, typename Output, typename Operator, typename ...Args>
  inline void gauss(std::vector<IntegrationPoint> const &ips, Output &out, Operator && op, Args... args)
  {
    for (auto const &ip : ips)
    {
      out += ip.weight*ip.jdet*op(ip,args...);
    }
  }
}