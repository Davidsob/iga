
#include "splines/Algorithms.h"
#include "splines/BSpline.h"
#include "splines/utils/VectorOperations.h"

#include <iostream>
#include <vector>
#include <string>



void curveTest()
{
  using namespace vector_ops;
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  BSplineCurve curve;

  int p = 3;
  typename decltype(curve)::vector knot{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
  typename decltype(curve)::matrix cpts{
    {0.0, 0.0},
    {1.0, 1.0},
    {2.0, 1.2},
    {3.0, 0.0}
  };

  curve.p = p;
  curve.knot = knot;
  curve.Q = cpts;
  std::cout << "point on curve = " << spline_ops::CurvePoint(0.5, curve) << std::endl;
  std::cout << "derivs = " << spline_ops::CurveDerivatives(0.5, 2, curve) << std::endl;
}

void surfaceTest()
{
  using namespace vector_ops;
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  BSplineSurface surf;

  int p = 3; int q = 1;
  typename decltype(surf)::vector uknot{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
  typename decltype(surf)::vector vknot{0.0, 0.0, 1.0, 1.0};
  // typename decltype(surf)::matrix cpts{
  //   {0.0, 0.0, 0.0},
  //   {1.0, 1.0, 0.0},
  //   {2.0, 1.2, 0.0},
  //   {3.0, 0.0, 0.0}, 

  //   {0.0, 0.0, 1.0},
  //   {1.0, 1.0, 1.0},
  //   {2.0, 1.2, 1.0},
  //   {3.0, 0.0, 1.0}
  // };

  using matrix = typename decltype(surf)::matrix;
  std::vector<matrix> cpts {
    {{0.0, 0.0, 0.0},
    {1.0, 1.0, 0.0},
    {2.0, 1.2, 0.0},
    {3.0, 0.0, 0.0}}, 

    {{0.0, 0.0, 1.0},
    {1.0, 1.0, 1.0},
    {2.0, 1.2, 1.0},
    {3.0, 0.0, 1.0}}
  };

  surf.p = p;
  surf.q = q;
  surf.uknot = uknot;
  surf.vknot = vknot;
  surf.Q = cpts;
  std::cout << "point on surface= " << spline_ops::SurfacePoint(0.5, 0.5, surf) << std::endl;
  auto ders = spline_ops::SurfaceDerivatives(0.5, 0.5, 1, surf);
  auto Su = ders[1][0];
  auto Sv = ders[0][1];
  auto S  = ders[0][0];

  std::cout << "S = " << S << std::endl;
  std::cout << "Su = " << Su << std::endl;
  std::cout << "Sv = " << Sv << std::endl;
}

int main(int argc, char **argv)
{
  std::cout << "*** B-Spline Main ***" << std::endl;
  curveTest();
  surfaceTest();
  return 0;
}
