#pragma once

#include <Eigen/Dense>

namespace quadrature
{

  template<typename T>
  T zero();

  template<>
  Eigen::MatrixXd zero() { return Eigen::MatrixXd::Zero(); }

  template<typename IntegrationPoint, typename Output, typename Operator>
  void gauss(std::vector<IntegrationPoint> const &ips, Output &out, Operator && op)
  {
    out = zero<Output>();
    for (auto const &ip : ips)
      sum += op(ip)*ip.w*ip.jac;
    return sum;
  }
}