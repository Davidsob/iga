#pragma once

#include <Eigen/Dense>

#include <vector>

namespace convert
{
  template<typename From, typename To>
  To _converter(From const &in);

  template<typename To, typename From>
  To to(From const &in) { return _converter<From,To>(in);}

  template<>
  Eigen::VectorXd _converter(std::vector<double> const &in)
  {
    return Eigen::Map<const Eigen::VectorXd>(&in[0],in.size());
  }

  template<>
  Eigen::VectorXd const _converter(std::vector<double> const &in)
  {
    return Eigen::Map<const Eigen::VectorXd>(&in[0],in.size());
  }

  template<>
  Eigen::RowVectorXd const _converter(std::vector<double> const &in)
  {
    return Eigen::Map<const Eigen::RowVectorXd>(&in[0],in.size());
  }

  template<>
  Eigen::RowVectorXd _converter(std::vector<double> const &in)
  {
    return Eigen::Map<const Eigen::RowVectorXd>(&in[0],in.size());
  }

  template<>
  Eigen::MatrixXd _converter(std::vector<std::vector<double>> const &in)
  {
    Eigen::MatrixXd out(in.size(),in[0].size());
    size_t row{0};
    for (auto const &x : in) out.row(row++) = to<Eigen::VectorXd>(x);
    return out;
  }

  // template<>
  // Eigen::VectorXd const _converter(std::vector<double> const &in)
  // {
  //   return Eigen::Map<const Eigen::VectorXd>(&in[0],in.size());
  // }

}