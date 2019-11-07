#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "utils/VectorOperations.h"

#ifndef FLT_TOL
#define FLT_TOL 1e-8
#endif

using namespace vector_ops;


namespace algo
{
  template<typename T>
  inline bool equal(T const &a, T const &b)
  { 
    return std::abs(a-b) <= FLT_TOL;
  }

  inline void normalizeKnot(std::vector<double> &knot)
  {
    auto min = knot[0];
    auto max = knot.back();
    auto L = max-min; 
    std::for_each(knot.begin(), knot.end(), [min,L](auto &x) { x = (x-min)/L; });
  }

  /**
   * @brief      { finds the interval in the knot vector which contains the given parametric coordinate }
   *
   * @param[in]  u     { parametric coordinate }
   * @param[in]  p     { polynomial order of spline }
   * @param      knot  The knot
   *
   * @return     { index of interval containing parametric coordinate.
   *               If the coordinate is not contained in the know vector then returns size_t(-1)
   *              }
   */
  inline size_t FindSpan(double u, size_t p, std::vector<double> const &knot)
  {
    auto m = knot.size()-1;
    auto n = m-p-1; 

    if (equal(u,knot[n+1])) return n; // special case

    if (u < knot[0] || u > knot.back()) {
      return size_t(-1); // special case
    }
    
    // binary search 
    auto low = p; auto high = n+1; auto mid = (low+high)/2;
    
    while (u < knot[mid] || u >= knot[mid+1])
    {
      if (u < knot[mid]) high = mid;
      else low = mid;
      mid = (low + high)/2;
    }

    return mid;
  }
 
  /**
   * @brief      Computes non-zero b-spline basis functions given coordinate in the knot span
   *
   * @param[in]  u     { span coordinate }
   * @param[in]  i     { level of shape function }
   * @param[in]  p     { polynomial order }
   * @param      knot  The knot
   *
   * @return     { non-zero basis functions}
   */
  std::vector<double>
  inline BasisFunctions(double u, int i, int p, std::vector<double> const &knot)
  {
    std::vector<double> N(p+1,0);
    if (i < p) return N;

    double saved{0.0};
    double tmp{0.0};
    std::vector<double> left(N);
    std::vector<double> right(N);

    N[0] = 1.0; 
    for (int j = 1; j <= p ; j++)
    {
      left[j] = u-knot[i+1-j]; 
      right[j] = knot[i+j]-u; 
      saved = 0.0; 
      for (int r = 0; r < j; r++)
      {
        tmp = N[r]/(right[r+1] + left[j-r]);
        N[r] = saved + right[r+1]*tmp;
        saved = left[j-r]*tmp;
      }
      N[j] = saved;
    }
    return N;
  }

  std::vector<std::vector<double>>
  inline BasisFunctionDerivatives(double u, int i, int p, std::vector<double> const &knot)
  {
    using matrix = std::vector<std::vector<double>>;

    int m = knot.size()-1;
    int n = m-p-1;

    double saved{0.0};
    double tmp{0.0};
    std::vector<double> left(p+1,0.0);
    std::vector<double> right(left);
    matrix ndu(n+1, std::vector<double>(left));
    matrix ders(ndu);
    matrix a(ndu);

    ndu[0][0] = 1.0;
    for (int j = 1; j <= p ; j++)
    {
      left[j] = u-knot[i+1-j]; 
      right[j] = knot[i+j]-u; 
      saved = 0.0; 
      for (int r = 0; r < j; r++)
      {
        ndu[j][r] = right[r+1]+left[j-r]; // lower triangle
        tmp = ndu[r][j-1]/ndu[j][r];

        ndu[r][j] = saved + right[r+1]*tmp; // upper triangle
        saved = left[j-r]*tmp;
      }
      ndu[j][j] = saved;
    }

    for (int j = 0; j <= p; j++)
      ders[0][j] = ndu[j][p];

    // compute derivatives as in eq. 2.9 of the nurbs book
    double d{0.0};
    int s1{0}; int s2{1};
    for (int r = 0; r <= p; r++)
    {
      s1 = 0; s2 = 1;
      a[0][0] = 1.0;
      // loop to compute kth der
      for (int k = 1; k <= n; k++)
      {
        d = 0.0;
        auto rk = r-k; auto pk = p-k;
        if (r >= k)
        {
          a[s2][0] = a[s1][0]/ndu[pk+1][rk];
          d = a[s2][0]*ndu[rk][pk];
        }
 
        auto j1 = (rk >= -1) ? 1 : -rk;
        auto j2 = (r-1 <= pk) ? k-1 : p-r;

        for (int j = j1; j <= j2; j++)
        {
          a[s2][j] = (a[s1][j] - a[s1][j-1])/ndu[pk+1][rk+j];
          d += a[s2][j]*ndu[rk+j][pk];
        }

        if (r <= pk)
        {
          a[s2][k] = -a[s1][k-1]/ndu[pk+1][r];
          d += a[s2][k]*ndu[r][pk];
        }

        ders[k][r] = d;
        std::swap(s1, s2); // switch rows
      }
    }
    // multiply by correction factors
    int r = p;
    for (int k = 1; k <= n; k++)
    {
      for (int j = 0; j <= p; j++)
      {
        ders[k][j] *= r;
      }
      r *= (p-k);
    }

    return ders;
  }

  double
  inline BasisFunction(double u, int i, int p, std::vector<double> const &knot)
  {
    int m = knot.size()-1;
    // compute the basis function Nip
    if ((i == 0 && equal(u,knot[0])) ||
        (i == m-p-1 && equal(u,knot[m])))
    {
      return 1.0;
    }

    if (u < knot[i] || u >= knot[i+p+1])
    {
      return 0.0; // local property
    }

    std::vector<double> N(p+1,0);
    for (int j = 0; j <= p; j++)
    {
      if (u >= knot[i+j] && u < knot[i+j+1]) 
      {
        N[j] = 1.0;
      }
    }

    double saved{0.0};
    for (int k = 1; k <= p; k++)
    {
      if (equal(N[0], 0.0))
      {
        saved = 0.0;
      } else {
        saved = ((u-knot[i])*N[0])/(knot[i+k]-knot[i]);
      }

      for (int j = 0; j < p-k+1; j++)
      {
        auto left = knot[i+j+1];
        auto right = knot[i+j+k+1];
        if (equal(N[j+1],0.0))
        {
          N[j] = saved;
          saved = 0.0;
        } else {
          auto tmp = N[j+1]/(right-left);
          N[j] = saved + (right-u)*tmp;
          saved = (u-left)*tmp;
        }
      }
    }
    return N[0];
  }

  std::vector<double> 
  inline BasisFunctionDerivative(double u, int i, int p, std::vector<double> const &knot)
  {
    using matrix = std::vector<std::vector<double>>;

    // int m = knot.size()-1;
    // int n = m-p-1;
    int n = p;

    std::vector<double> ders(n+1,0);
    if (u < knot[i] || u >= knot[i+p+1]) // local property
    {
      return ders;
    }

    matrix N(p+1,std::vector<double>(p+1,0.0));
    for (int j = 0; j <= p; j++)
    {
      if (u >= knot[i+j] && u < knot[i+j+1])
      {
        N[j][0] = 1.0;
      }
    }

    // compute triangular table
    double saved{0.0};
    for (int k = 1; k <= p; k++)
    {
      if (equal(N[0][k-1], 0.0)) 
      {
        saved = 0.0;
      } else {
        saved = ((u-knot[i])*N[0][k-1])/(knot[i+k]-knot[i]);
      }

      for (int j = 0; j < p-k+1; j++)
      {
        auto left = knot[i+j+1];
        auto right = knot[i+j+k+1];
        if (equal(N[j+1][k-1],0.0))
        {
          N[j][k] = saved;
          saved = 0.0; 
        } else {
          auto tmp = N[j+1][k-1]/(right-left);
          N[j][k] = saved + (right-u)*tmp;
          saved = (u-left)*tmp;
        }
      }
    }

    std::vector<double> Nd(ders.size(), 0.0);
    ders[0] = N[0][p]; // function value
    for (int k = 1; k <= n; k++)
    {
      for (int j = 0; j <= k; j++)
      {
        Nd[j] = N[j][p-k]; // load column
      }

      for (int jj = 1; jj <= k; jj++) // compute table width
      {
        if (equal(Nd[0],0.0))
        {
          saved = 0.0;
        } else {
          saved = Nd[0]/(knot[i+p-k+jj] - knot[i]);
        }
        for (int j = 0; j < k-jj+1; j++)
        {
          auto left = knot[i+j+1];
          auto right = knot[i+j+p+jj+1];
          if (equal(Nd[j+1],0.0))
          {
            Nd[j] = (p-k+jj)*saved;
            saved = 0.0;
          } else {
            auto tmp = Nd[j+1]/(right-left);
            Nd[j] = (p-k+jj)*(saved -tmp);
            saved = tmp;
          }
        }
      }
      ders[k] = Nd[0]; // the kth derivative
    }

    return ders;
  }
}