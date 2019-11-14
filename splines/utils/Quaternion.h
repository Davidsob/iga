#pragma once

#include "VectorOperations.h"

#include <iostream>
#include <vector>
#include <cmath>

using namespace vector_ops;

class Quaternion
{
public:
  Quaternion(double const &x,std::vector<double> const &u)
    : s(x), v(u)
  {}

  Quaternion(Quaternion const &copy)
  {
    s = copy.s;
    v = copy.v;
  }

  Quaternion(Quaternion &&tmp)
  {
    s = tmp.s;
    v = std::move(tmp.v);
  }

  virtual ~Quaternion() = default;

  void setScalarComponent(double const &x) { s = x; }
  double const &scalarComponent() const { return s; }

  void setVectorComponent(std::vector<double> const &x) { v = x; }
  std::vector<double> const &vectorComponent() const { return v; }

  friend std::ostream& operator<<(std::ostream &os, Quaternion const &x);

protected:
  double s; // scalar part
  std::vector<double> v; // vector part
};

inline std::ostream &operator<<(std::ostream &os, Quaternion const &x)
{
  os << "< " << x.scalarComponent() << " ";
  for(auto e : x.vectorComponent()) os << e << " ";
  os << ">";
  return os; 
}

inline Quaternion operator+(Quaternion const &a, Quaternion const &b)
{
  return Quaternion(a.scalarComponent() + b.scalarComponent(),
                    a.vectorComponent() + b.vectorComponent());
}

inline void operator+=(Quaternion &a, Quaternion const &b)
{
  a.setScalarComponent(a.scalarComponent() + b.scalarComponent());
  a.setVectorComponent(a.vectorComponent() + b.vectorComponent());
}

inline Quaternion operator-(Quaternion const &a, Quaternion const &b)
{
  return Quaternion(a.scalarComponent() - b.scalarComponent(),
                    a.vectorComponent() - b.vectorComponent());
}

inline void operator-=(Quaternion &a, Quaternion const &b)
{
  a.setScalarComponent(a.scalarComponent() - b.scalarComponent());
  a.setVectorComponent(a.vectorComponent() - b.vectorComponent());
}

template<typename T>
inline Quaternion operator*(Quaternion const &a, T const &c)
{
  return Quaternion(c*a.scalarComponent(), c*a.vectorComponent());
}

template<typename T>
inline Quaternion operator*(T const &c, Quaternion const &a)
{
  return a*c;
}

template<typename T>
inline void operator*=(Quaternion &a, T const &c)
{
  a.setScalarComponent(a.scalarComponent()*c);
  a.setVectorComponent(a.vectorComponent()*c);
}

inline Quaternion operator*(Quaternion const &a, Quaternion const &b)
{
  auto const &s1 = a.scalarComponent();
  auto const &s2 = b.scalarComponent();
  auto const &v1 = a.vectorComponent();
  auto const &v2 = b.vectorComponent();

  return Quaternion(
            s1*s2 - dot(v1,v2),
            s1*v2 + s2*v1 + cross(v1,v2)
          );
}

template<typename T>
inline Quaternion operator/(Quaternion const &a, T const &c)
{
  return Quaternion(a.scalarComponent()/c, a.vectorComponent()/c);
}

template<typename T>
inline void operator/=(Quaternion &a, T const &c)
{
  a.setScalarComponent(a.scalarComponent()/c);
  a.setVectorComponent(a.vectorComponent()/c);
}

inline Quaternion conj(Quaternion const &a)
{
  return Quaternion(a.scalarComponent(), -a.vectorComponent());
}

inline double norm(Quaternion const &a)
{
  return std::sqrt(Quaternion(conj(a)*a).scalarComponent());
}

inline void normalize(Quaternion &a)
{
  a /= norm(a);
}

inline Quaternion inverse(Quaternion const &a)
{
  return conj(a)/std::pow(norm(a),1);
}

inline Quaternion rotate(Quaternion const &q, Quaternion const &p)
{
  return q*(p*inverse(q));
}
