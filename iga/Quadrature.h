#pragma once

#include <Eigen/Dense>

#include <type_traits>

namespace quadrature
{

  template<typename T>
  inline T zero(T const &);

  template<>
  inline Eigen::MatrixXd zero(Eigen::MatrixXd const &mat) { return Eigen::MatrixXd::Zero(mat.rows(),mat.cols()); }

  template<>
  inline Eigen::VectorXd zero(Eigen::VectorXd const &vec) { return Eigen::VectorXd::Zero(vec.size()); }

  template<typename IntegrationPoint, typename Output, typename Operator>
  inline void gauss(std::vector<IntegrationPoint> const &ips, Output &out, Operator && op)
  {
    Output sum = zero(out);
    for (auto const &ip : ips)
    {
      sum += ip.weight*ip.jdet*op(ip);
    }
    out += sum;
  }
}