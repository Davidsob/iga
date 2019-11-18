#include "SplineModifiers.h"

#include "Algorithms.h"
#include "BSpline.h"
#include "Nurbs.h"
#include "utils/VectorOperations.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace spline_ops
{
  void elevate(size_t t, BSplineCurve &curve)
  {
    auto const p = curve.p;
    auto const &Pw = curve.Q;
    auto const &U = curve.knot;
    auto const n = U.size()-p-1;

    auto const m = n+p+1;
    auto const ph = p+t;
    auto ph2 = ph/2;

    // compute bezier degree elevation coefficeints 
    std::vector<vector<double>> bezalfs(p+t+1, std::vector<double>(p+1,0));
    bezalfs[0][0] = bezalfs[ph][p] = 1.0;
    for (size_t i = 1; i <= ph2 ; i++)
    {
        auto inv = 1.0/algo::
    }
  }

  void elevate(size_t t, NurbsCurve &curve)
  {

  }
}
