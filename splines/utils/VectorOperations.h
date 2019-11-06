#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <ostream>
#include <fstream>
#include <vector>

namespace vector_ops 
{
  template<typename T>
  std::ostream &operator<<(std::ostream &os, std::vector<T> const &x)
  {
    os << "< ";
    for(auto e : x) os << e << " ";
    os << ">";
    return os; 
  }

  template<typename T>
  std::vector<T> operator+(std::vector<T> const &a, std::vector<T> const &b)
  {
    static auto plus = [](T const &x, T const &y) { return x+y; };
    std::vector<T> c(a.size(),T());
    std::transform(a.begin(), a.end(), b.begin(),c.begin(), plus);
    return c;
  }

  template<typename T>
  void operator+=(std::vector<T> &a, std::vector<T> const &b)
  {
    static auto plus = [](T const &x, T const &y) { return x+y; };
    std::transform(a.begin(), a.end(), b.begin(),a.begin(), plus);
  }

  template<typename T>
  std::vector<T> operator-(std::vector<T> const &a, std::vector<T> const &b)
  {
    static auto minus = [](T const &x, T const &y) { return x-y; };
    std::vector<T> c(a.size(),T());
    std::transform(a.begin(), a.end(), b.begin(),c.begin(), minus);
    return c;
  }

  template<typename T>
  std::vector<T> operator-(std::vector<T> const &a)
  {
    static auto minus = [](T const &x) { return -x; };
    std::vector<T> c(a.size(),T());
    std::transform(a.begin(), a.end(), c.begin(), minus);
    return c;
  }

  template<typename T>
  void operator-=(std::vector<T> &a, std::vector<T> const &b)
  {
    static auto minus = [](T const &x, T const &y) { return x-y; };
    std::transform(a.begin(), a.end(), b.begin(),a.begin(),minus);
  }

  template<typename T>
  std::vector<T> operator*(std::vector<T> const &x, std::vector<T> const &y)
  {
    std::vector<T> out(x);
    std::transform(out.begin(), out.end(), y.begin(),out.begin(),std::multiplies<T>());
    return out;
  }

  template<typename T, typename U>
  std::vector<T> operator*(U const &c, std::vector<T> const &x)
  {
    std::vector<T> y(x.size(),T());
    std::transform(x.begin(), x.end(), y.begin(), [&c](T const xi){return c*xi;});
    return y;
  }

  template<typename T, typename U>
  std::vector<T> operator*(std::vector<T> const &x, U const &c)
  {
    return c*x;
  }

  template<typename T, typename U>
  void operator*=(std::vector<T> &x, U const &c)
  {
    std::transform(x.begin(), x.end(), x.begin(), [&c](T const xi){return c*xi;});
  }

  template<typename T, typename U>
  std::vector<T> operator/(std::vector<T> const &x, U const &c)
  {
    std::vector<T> y(x.size(),T());
    std::transform(x.begin(), x.end(), y.begin(), [&c](T const xi){return xi/c;});
    return y;
  }

  template<typename T, typename U>
  void operator/=(std::vector<T> &x, U const &c)
  {
    std::transform(x.begin(), x.end(), x.begin(), [&c](T const xi){return xi/c;});
  }

  template<typename T>
  T dot(std::vector<T> const &a, std::vector<T> const &b)
  {
    std::vector<T> c(a.size(),T());
    std::transform(a.begin(), a.end(), b.begin(),c.begin(),std::multiplies<T>());
    T zero(0);
    T sum = std::accumulate(c.begin(), c.end(), zero);
    return sum;
  }

  template<typename T>
  std::vector<T> dot(std::vector<T> const &a, std::vector<std::vector<T>> const &b)
  {
    std::vector<std::vector<T>> c(a.size(), std::vector<T>(b[0].size(),0));
    std::transform(a.begin(), a.end(), b.begin(),c.begin(),
      [](T const &a, std::vector<T> const &b)
      {
        return a*b;
      }
    );

    std::vector<T> zero(b[0].size(),0);
    auto sum = std::accumulate(c.begin(), c.end(), zero,
      [](auto const &a, auto const &b)
      {
        return a+b;
      }
    );
    return sum ;
  }

  template<typename T>
  std::vector<T> cross(std::vector<T> const &a, std::vector<T> const &b)
  {
    std::vector<T> c(a.size(),T());
    c[0] = a[1]*b[2]-a[2]*b[1];
    c[1] = a[2]*b[0]-a[0]*b[2];
    c[2] = a[0]*b[1]-a[1]*b[0];
    return c;
  }

  template<typename T>
  T norm(std::vector<T> const &a)
  {
    return sqrt(dot(a,a));
  }

  template<typename T>
  T inf_norm(std::vector<T> const &a)
  {
    T mx = std::abs(a[0]);
    for (size_t i = 0; i < a.size(); i++)
    {
      mx = std::max(mx, std::abs(a[i]));
    }
    return mx;
  }

  template<typename T>
  std::vector<T> normalize(std::vector<T> const &x)
  {
    auto n = norm(x);
    auto fact = ( n < 1e-8) ? T(1.0) : T(1.0)/n;
    return fact*x;
  }
}
