#ifndef GeometryUtils_h
#define GeometryUtils_h

#include "MatrixTypes.h"

#include <iostream>

namespace GeometryUtils
{
  
  struct Plane
  {
    explicit Plane(StaticVectorR<2> const &x, StaticVectorR<2> const &n)
      : x(x), n(n)
    { }
    
    virtual ~Plane() = default;

    StaticVectorR<2> x, n;  
  };

  struct Circle
  {
    explicit Circle(StaticVectorR<2> const &x, double R = 1e-8)
      : x(x), r(R)
    { }
    
    virtual ~Circle() = default;

    StaticVectorR<2> x;
    double r;
  };

  template<typename T>
  bool on_plane(T const &x, Plane const &plane, double tol=1e-8)
  {
    return std::abs(plane.n.dot(x - plane.x)) <= tol;
  }

  template<typename T>
  bool in_halfspace(T const &x, Plane const &plane, double tol=1e-8)
  {
    return plane.n.dot(x - plane.x) >= 0.0 || on_plane(x, plane, tol);
  }

  template<typename T>
  bool edge_on_plane(T const &e, Plane const &plane, double tol=1e-8)
  {
    bool on = true;
    for (auto &x : e.x)
    {
      on &= on_plane(x, plane, tol);
    }
    on &= (std::abs(e.normal.dot(plane.n) - 1.0) <= tol);
    return on;
  }

  StaticVectorR<2>
  normal2(StaticVectorR<2> const &a, StaticVectorR<2> const &b);

  template<typename T>
  bool in_circle(T const &x, Circle const &circle, double tol=1e-8)
  {
    return (x - circle.x).norm() <= circle.r;
  }
  
}

#endif // GeometryUtils_h
