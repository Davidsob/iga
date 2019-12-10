#pragma once

#include "MatrixTypes.h"

// traits
template<typename T>
struct is_matrix_type
{
  static constexpr bool value() { return false; }
};

template<>
struct is_matrix_type<double>
{
  static constexpr bool value() { return false; }
};

template<>
struct is_matrix_type<DynamicMatrixR>
{
  static constexpr bool value() { return true; }
};

template<>
struct is_matrix_type<SparseMatrixR>
{
  static constexpr bool value() { return true; }
};

template<>
struct is_matrix_type<StaticMatrixR<1,1>>
{
  static constexpr bool value() { return true; }
};

template<>
struct is_matrix_type<StaticMatrixR<2,2>>
{
  static constexpr bool value() { return true; }
};

template<>
struct is_matrix_type<StaticMatrixR<4,4>>
{
  static constexpr bool value() { return true; }
};

template<>
struct is_matrix_type<StaticMatrixR<8,8>>
{
  static constexpr bool value() { return true; }
};

template<>
struct is_matrix_type<StaticMatrixR<16,16>>
{
  static constexpr bool value() { return true; }
};

template<typename T>
struct is_vector_type
{
  static constexpr bool value() { return !is_matrix_type<T>::value(); }
};

template<>
struct is_vector_type<StaticVectorR<1>>
{
  static constexpr bool value() { return true; }
};

template<>
struct is_vector_type<double>
{
  static constexpr bool value() { return false; }
};
