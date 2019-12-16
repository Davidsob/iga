#pragma once

#include <Eigen/Dense>

#include <vector>

#include "utils/MatrixTypes.h"

namespace convert
{
  template<typename From, typename To, typename ...Args>
  inline To _converter(From const &in, Args &...args);

  template<typename To, typename From, typename ...Args>
  inline To to(From const &in, Args &...args) { return _converter<From,To>(in,args...);}

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

  template<>
  inline std::vector<std::vector<double>> _converter(Eigen::MatrixXd const &in)
  {
    std::vector<std::vector<double>> out;
    for (int i = 0; i < in.rows(); i++) out.push_back(to<std::vector<double>>(Eigen::RowVectorXd(in.row(i))));
    return out;
  }

  template<>
  inline std::vector<double> _converter(Eigen::MatrixXd const &in)
  {
    std::vector<double> out;
    for (int i = 0; i < in.rows(); i++) out.push_back(in(i,0));
    return out;
  }

  template<>
  inline std::vector<Eigen::RowVectorXd> _converter(Eigen::MatrixXd const &in)
  {
    std::vector<Eigen::RowVectorXd> out;
    for (int i = 0; i < in.rows(); i++) out.push_back(Eigen::RowVectorXd(in.row(i)));
    return out;
  }

  template<>
  inline std::vector<Eigen::RowVector3d> _converter(Eigen::MatrixXd const &in)
  {
    std::vector<Eigen::RowVector3d> out;
    for (int i = 0; i < in.rows(); i++) out.push_back(Eigen::RowVector3d(in.row(i)));
    return out;
  }

  template<>
  inline std::vector<StaticRowVectorR<6>> _converter(Eigen::MatrixXd const &in)
  {
    std::vector<StaticRowVectorR<6>> out;
    for (int i = 0; i < in.rows(); i++) out.push_back(StaticRowVectorR<6>(in.row(i)));
    return out;
  }

  template<>
  inline void _converter(Eigen::VectorXd const &in,
                         std::vector<size_t> const &dof,
                         std::vector<Triplet> &out)
  {
    for (int i = 0; i < in.size(); i++)
      out.push_back(Triplet(dof[i], 0, in[i]));
  }

  template<>
  inline void _converter(Eigen::MatrixXd const &in,
                         std::vector<size_t> const &rdof,
                         std::vector<size_t> const &cdof,
                         std::vector<Triplet> &out)
  {
    SparseMatrixR const &sp = in.sparseView();
    for (int k = 0; k < sp.outerSize(); k++)
    {
      for (SparseMatrixR::InnerIterator it(sp,k); it; ++it)
      {
        out.push_back(Triplet(rdof[it.row()], cdof[it.col()], it.value()));
      }
    }
  }
}