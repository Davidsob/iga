#pragma once

#include <iostream>

class BSplineCurve;
class BSplineSurface;
class BSplineSolid;
class NurbsCurve;
class NurbsSurface;
class NurbsSolid;

namespace spline_ops
{
  // degree reduction 
  bool reduce(BSplineCurve &, double tolerance = 1e-3);
  void reduce(NurbsCurve   &, double tolerance = 1e-3);

  void reduce(size_t direction, BSplineSurface &, double tolerance = 1e-3);
  void reduce(size_t direction, NurbsSurface   &, double tolerance = 1e-3);
}