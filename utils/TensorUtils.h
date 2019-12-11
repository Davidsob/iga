#ifndef TensorUtils_h
#define TensorUtils_h

#include "MatrixTypes.h"

struct TensorUtils
{
  static const StaticMatrixR<3,3> eye;
  static const StaticMatrixR<6,6> eye6;
  static const StaticMatrixR<9,9> eye9;
  static const size_t ik[6];
  static const size_t jl[6];

  static
  StaticVectorR<4>
  stressToVoigt(StaticMatrixR<3,3> const &x)
  {
    StaticVectorR<4> y; y << x(0,0), x(1,1), x(2,2), x(0,1);
    return y;
  }

  static
  StaticVectorR<4>
  strainToVoigt(StaticMatrixR<3,3> const &x)
  {
    StaticVectorR<4> y; y << x(0,0), x(1,1), x(2,2), 2.0*x(0,1);
    return y;
  }

  static
  StaticMatrixR<3,3>
  stressFromVoigt(StaticVectorR<4> const &x)
  {
    StaticMatrixR<3,3> y; y.setZero();
    y(0,0) = x[0];
    y(1,1) = x[1];
    y(2,2) = x[2];
    y(0,1) = x[3];
    y(1,0) = x[3];
    return y;
  }

  static
  StaticMatrixR<3,3>
  strainFromVoigt(StaticVectorR<4> const &x)
  {
    StaticMatrixR<3,3> y; y.setZero();
    y(0,0) = x[0];
    y(1,1) = x[1];
    y(2,2) = x[2];
    y(0,1) = 0.5*x[3];
    y(1,0) = 0.5*x[3];
    return y;
  }

};

#endif // TensorUtils_h
