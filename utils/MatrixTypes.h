#ifndef MatrixTypes_h
#define MatrixTypes_h

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iostream>

template<int N, int M>
using StaticMatrixR = Eigen::Matrix<double, N, M>;

template<int N>
using StaticVectorR = Eigen::Matrix<double, N, 1>;

template<int N>
using StaticRowVectorR = Eigen::Matrix<double, 1, N>;

using DynamicMatrixR = Eigen::MatrixXd;
using DynamicVectorR = Eigen::VectorXd;
using DynamicRowVectorR = Eigen::RowVectorXd;

struct Triplet : public Eigen::Triplet<double>
{
  explicit Triplet(size_t r, size_t c, double x)
    : Eigen::Triplet<double>(r,c,x) {}

  void setRow(size_t i) { m_row = i; }
  void setCol(size_t i) { m_col = i; }
  void setValue(double x) { m_value = x; }
  void transpose() { std::swap(m_row, m_col); }

  friend std::ostream & operator<<(std::ostream &os, Triplet const &t)
  {
    os << "{" << t.row() << ", " << t.col() << ", " << t.value() << "}";
    return os;
  }
};

using SparseMatrixR = Eigen::SparseMatrix<double>;

#endif //MatrixTypes_h
