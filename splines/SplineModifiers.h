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
  // knot insertion algorithms
  void insertKnot(double u, int r, BSplineCurve &);
  void insertKnot(double uv, int r, int direction, BSplineSurface &);
  void insertKnot(double uvw, int r, int direction, BSplineSolid &);

  void insertKnot(double u, size_t rep, NurbsCurve &);
  void insertKnot(double uv, size_t rep, size_t direction, NurbsSurface &);
  void insertKnot(double uvw, size_t rep, size_t direction, NurbsSolid &s);

  // knot refinement routines
  template<typename Shape>
  void midpointRefinement(size_t level, Shape&); // uniform refinement - all directions

  template<typename Shape>
  void midpointRefinement(size_t level, size_t direction, Shape &);


  // degree elevation
  void elevate(size_t rep, BSplineCurve &);
  void elevate(size_t rep, NurbsCurve &);
}