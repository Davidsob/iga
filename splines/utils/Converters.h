#pragma once

#include <Eigen/Dense>

#include <vector>

namespace convert
{
  template<typename From, typename To>
  inline To _converter(From const &in);

  template<typename To, typename From>
  inline To to(From const &in) { return _converter<From,To>(in);}

  template<>
  inline Eigen::VectorXd _converter(std::vector<double> const &in)
  {
    return Eigen::Map<const Eigen::VectorXd>(&in[0],in.size());
  }

  template<>
  inline Eigen::VectorXd const _converter(std::vector<double> const &in)
  {
    return Eigen::Map<const Eigen::VectorXd>(&in[0],in.size());
  }

  template<>
  inline Eigen::RowVectorXd const _converter(std::vector<double> const &in)
  {
    return Eigen::Map<const Eigen::RowVectorXd>(&in[0],in.size());
  }

  template<>
  inline Eigen::RowVectorXd _converter(std::vector<double> const &in)
  {
    return Eigen::Map<const Eigen::RowVectorXd>(&in[0],in.size());
  }

  template<>
  inline std::vector<double> _converter(Eigen::VectorXd const &in)
  {
    return std::vector<double>(in.data(), in.data()+in.size());
  }

  template<>
  inline std::vector<double> const _converter(Eigen::VectorXd const &in)
  {
    return std::vector<double>(in.data(), in.data()+in.size());
  }

  template<>
  inline std::vector<double> _converter(Eigen::RowVectorXd const &in)
  {
    return std::vector<double>(in.data(), in.data()+in.size());
  }

  template<>
  inline std::vector<double> const _converter(Eigen::RowVectorXd const &in)
  {
    return std::vector<double>(in.data(), in.data()+in.size());
  }

  template<>
  inline Eigen::MatrixXd _converter(std::vector<std::vector<double>> const &in)
  {
    Eigen::MatrixXd out(in.size(),in[0].size());
    size_t row{0};
    for (auto const &x : in) out.row(row++) = to<Eigen::VectorXd>(x);
    return out;
  }
}