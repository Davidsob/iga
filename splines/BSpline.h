#pragma once

#include "utils/VectorOperations.h"

#include "Algorithms.h"

#include <iostream>
#include <fstream>

#include <string>


struct BSplineCurve
{
  using vector = std::vector<double>;
  using matrix = std::vector<std::vector<double>>; // need a point type

  size_t dim() const { return Q.empty() ? 0 : Q[0].size(); }
  int p; //polynomial order
  vector knot;
  matrix Q; 
};

struct BSplineSurface
{
  using vector = std::vector<double>;
  using matrix = std::vector<std::vector<double>>; // need a point type

  int p, q; //polynomial order
  vector uknot, vknot; // knot vectors
  std::vector<matrix> Q; // cpts in vector form {c0j, c1j, cij...} 
};

namespace spline_ops
{
  std::vector<double>
  CurvePoint(double u, BSplineCurve const &curve)
  {
    using namespace vector_ops;
    using point = typename BSplineCurve::matrix::value_type;

    auto span = algo::FindSpan(u, curve.p, curve.knot);
    std::vector<double> N = algo::BasisFunctions(u,span,curve.p, curve.knot);

    point C(curve.dim(),0.0);
    for (int i = 0; i <= curve.p; i++)
    {
      C += N[i]*curve.Q[span-curve.p+i];
    }

    return C;
  }

  std::vector<std::vector<double>>
  CurveDerivatives(double u, int order, BSplineCurve const &curve)
  {
    using namespace vector_ops;
    using point = typename BSplineCurve::matrix::value_type;

    auto span = algo::FindSpan(u, curve.p, curve.knot);
    auto ders = algo::BasisFunctionDerivatives(u,span,curve.p,curve.knot);

    auto du = std::min(order,curve.p);
    decltype(ders) C(du+1, point(curve.dim(),0.0));
    for (int k = 0; k <= du; k++)
    {
      for (int j = 0; j <= curve.p; j++)
      {
        C[k] += ders[k][j]*curve.Q[span-curve.p+j];
      }
    }

    return C;
  }

  std::vector<double>
  SurfacePoint(double u, double v, BSplineSurface const &surf)
  {
    using namespace vector_ops;
    using point = typename BSplineSurface::matrix::value_type;

    auto uSpan = algo::FindSpan(u, surf.p, surf.uknot);
    auto vSpan = algo::FindSpan(v, surf.q, surf.vknot);
    std::vector<double> Nu = algo::BasisFunctions(u,uSpan,surf.p, surf.uknot);
    std::vector<double> Nv = algo::BasisFunctions(v,vSpan,surf.q, surf.vknot);

    point S(3,0.0);
    std::vector<point> tmp(surf.q+1,S);

    for (int i = 0; i <= surf.q; i++)
    {
      for (int j = 0; j <= surf.p; j++)
      {
        auto a = uSpan-surf.p+j;
        auto b = vSpan-surf.q+i;
        tmp[i] += Nu[j]*surf.Q[b][a];
      }
    }

    for (int i = 0; i <= surf.q; i++)
      S += Nv[i]*tmp[i];

    return S;
  }

  std::vector<std::vector<std::vector<double>>>
  SurfaceDerivatives(double u, double v, int order, BSplineSurface const &surf)
  {
    using namespace vector_ops;
    using point = typename BSplineSurface::matrix::value_type;
    using matrix = typename BSplineSurface::matrix;

    auto uspan = algo::FindSpan(u, surf.p, surf.uknot);
    auto vspan = algo::FindSpan(v, surf.q, surf.vknot);
    auto Nu = algo::BasisFunctionDerivatives(u,uspan,surf.p,surf.uknot);
    auto Nv = algo::BasisFunctionDerivatives(v,vspan,surf.q,surf.vknot);

    auto du = std::min(order,surf.p);
    auto dv = std::min(order,surf.q);

    point zero(3,0);
    auto dskl = std::max(du,dv);
    std::vector<matrix> Skl(dskl+1, matrix(dskl+1, zero));

    std::vector<point> tmp(surf.q+1,zero);
    for (int k = 0; k <= du; k++)
    {
      for (int s = 0; s <= surf.q; s++)
      {
        tmp[s] = zero;
        for (int r = 0; r <= surf.p; r++)
        {
          int a = uspan-surf.p+r;
          int b = vspan-surf.q+s;
          tmp[s] += Nu[k][r]*surf.Q[b][a];
        }
        auto dd = std::min(order-k,dv);
        for (int l = 0; l <= dd; l++)
        {
          Skl[k][l] = zero;
          for (int s = 0; s <= surf.q; s++)
            Skl[k][l] += Nv[l][s]*tmp[s];
        }
      }
    }
    return Skl;
  }

  template<typename T> void writeToFile;

  template<>
  void writeToFile(BSplineCurve const &c, std::string const &file, std::string const &dir="../output")
  {

  }
}