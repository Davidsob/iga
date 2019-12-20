#ifndef GeometryUtils_h
#define GeometryUtils_h

#include "MatrixTypes.h"

#include <iostream>

namespace GeometryUtils
{
  template<size_t Dim> 
  struct Plane
  {
    static constexpr size_t const dim = Dim;
    explicit Plane(StaticVectorR<Dim> const &x, StaticVectorR<Dim> const &n)
      : x(x), n(n)
    { }
    
    virtual ~Plane() = default;

    StaticVectorR<Dim> x, n;  
  };

  template<size_t Dim> 
  struct Circle
  {
    static constexpr size_t const dim = Dim;
    explicit Circle(StaticVectorR<Dim> const &x, double R = 1e-8)
      : x(x), r(R)
    { }
    
    virtual ~Circle() = default;

    StaticVectorR<Dim> x;
    double r;
  };

  template<typename T, typename Plane>
  bool on_plane(T const &x, Plane const &plane, double tol=1e-8)
  {
    return std::abs(plane.n.dot(x - plane.x)) <= tol;
  }

  template<typename T, typename Plane>
  bool in_halfspace(T const &x, Plane const &plane, double tol=1e-8)
  {
    return plane.n.dot(x - plane.x) >= 0.0 || on_plane(x, plane, tol);
  }

  template<typename T, typename Plane>
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

  template<typename T, typename Circle>
  bool in_circle(T const &x, Circle const &circle, double tol=1e-8)
  {
    return (x - circle.x).norm() <= circle.r;
  }
  
}

#endif // GeometryUtils_h
