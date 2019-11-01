#pragma once

#include "utils/VectorOperations.h"

#include "Algorithms.h"

#include <iostream>
#include <fstream>

#include <string>
#include <vector>

#include "BSpline.h"

struct NurbsCurve 
{
  using vector = std::vector<double>;
  using matrix = std::vector<std::vector<double>>; // need a point type

  size_t dim() const { return Q.empty() ? 0 : Q[0].size(); }
  int p; //polynomial order
  vector knot, weights;
  matrix Q; 
};

struct NurbsSurface 
{
  using vector = std::vector<double>;
  using matrix = std::vector<std::vector<double>>; // need a point type

  int qid(int iu, int iv) const { return iu + iv*(uknot.size()-p-1); }
  int p, q; //polynomial order
  vector uknot, vknot, weights; // knot vectors, weight vector
  matrix Q; // cpts in vector form {c0j, c1j, cij...} 
};

namespace spline_ops
{

  template<typename Nurb>
  void weightedControlPoints(Nurb const &c, typename Nurb::matrix &Qw)
  {
    using namespace vector_ops;
    // create weighted control points
    Qw = c.Q;
    for (size_t i = 0; i < c.weights.size(); i++)
    {
      double const &wi = c.weights[i];
      Qw[i] *= wi;
      Qw[i].push_back(wi);
    }
  }

  void weightedBSpline(NurbsCurve const &c, BSplineCurve &b)
  {
    b.p = c.p;
    b.knot = c.knot;
    weightedControlPoints(c,b.Q);
  }

  void weightedBSpline(NurbsSurface const &s, BSplineSurface &b)
  {
    b.p = s.p;
    b.q = s.p;
    b.uknot = s.uknot;
    b.vknot = s.vknot;
    weightedControlPoints(s,b.Q);
  }

  std::vector<double>
  CurvePoint(double u, NurbsCurve const &curve)
  {
    using namespace vector_ops;
    BSplineCurve b; weightedBSpline(curve, b);
    auto Cw = spline_ops::CurvePoint(u,b);
    Cw /= Cw.back(); 
    Cw.pop_back();
    return Cw;
  }

  inline int factorial(int k)
  {
    return k > 0 ? k*factorial(k-1) : 1.0;
  }

  inline double binomial(int n, int k)
  {
    return factorial(n)/(double(factorial(n-k)*factorial(k)));
  }

  std::vector<std::vector<double>>
  CurveDerivatives(double u, int order, NurbsCurve const &curve)
  {
    using namespace vector_ops;
    BSplineCurve b; weightedBSpline(curve, b);
    auto ders = spline_ops::CurveDerivatives(u,order,b);
    std::vector<std::vector<double>> dC;

    for (int k = 0; k <= order; k++)
    {
      auto v = ders[k];
      for (int i = 1; i <= k; i++)
      {
        v -= binomial(k,i)*ders[i].back()*dC[k-i];
      }

      v /= ders[0].back();
      v.pop_back();
      dC.push_back(v);
    }

    return dC;
  }

  std::vector<double>
  SurfacePoint(double u, double v, NurbsSurface const &surf)
  {
    using namespace vector_ops;
    BSplineSurface b; weightedBSpline(surf,b);
    std::cout << b.Q << std::endl;
    auto Sw = spline_ops::SurfacePoint(u,v,b);
    Sw /= Sw.back(); Sw.pop_back();
    return Sw;
  }

  // std::vector<std::vector<std::vector<double>>>
  // SurfaceDerivatives(double u, double v, int order, BSplineSurface const &surf)
  // {
  //   using namespace vector_ops;
  //   using point = typename BSplineSurface::matrix::value_type;
  //   using matrix = typename BSplineSurface::matrix;

  //   auto uspan = algo::FindSpan(u, surf.p, surf.uknot);
  //   auto vspan = algo::FindSpan(v, surf.q, surf.vknot);
  //   auto Nu = algo::BasisFunctionDerivatives(u,uspan,surf.p,surf.uknot);
  //   auto Nv = algo::BasisFunctionDerivatives(v,vspan,surf.q,surf.vknot);

  //   auto du = std::min(order,surf.p);
  //   auto dv = std::min(order,surf.q);

  //   point zero(3,0);
  //   auto dskl = std::max(du,dv);
  //   std::vector<matrix> Skl(dskl+1, matrix(dskl+1, zero));

  //   std::vector<point> tmp(surf.q+1,zero);
  //   for (int k = 0; k <= du; k++)
  //   {
  //     for (int s = 0; s <= surf.q; s++)
  //     {
  //       tmp[s] = zero;
  //       for (int r = 0; r <= surf.p; r++)
  //       {
  //         int a = uspan-surf.p+r;
  //         int b = vspan-surf.q+s;
  //         int idx = surf.qid(a,b);
  //         tmp[s] += Nu[k][r]*surf.Q[idx];
  //       }
  //       auto dd = std::min(order-k,dv);
  //       for (int l = 0; l <= dd; l++)
  //       {
  //         Skl[k][l] = zero;
  //         for (int s = 0; s <= surf.q; s++)
  //         {
  //           Skl[k][l] += Nv[l][s]*tmp[s];
  //         }
  //       }
  //     }
  //   }
  //   return Skl;
  // }

  void writeToFile(NurbsCurve const &c,std::string const &file_name, int level = 20)
  {
    using namespace vector_ops;

    // compute mesh
    std::vector<std::vector<int>> mesh;
    for (int i = 0; i < level; i++)
    {
      mesh.push_back({i,i+1});
    }
    // compute points on curve and write to file
    auto du = (c.knot.back() - c.knot[0])/level;
    auto u = 0.0;
    std::vector<std::vector<double>> pts;
    for (int i = 0; i <= level; i++)
    {
      pts.push_back(CurvePoint(u,c));
      u += du;
    }

    // write file
    std::ofstream file;
    file.open(file_name);
    file << mesh.size() << std::endl; // number of points in file
    for (auto &el : mesh) {
      for (auto &v : el) file << v << " ";
      file << std::endl;
    }

    file << pts.size() << std::endl;
    for (auto &p : pts) {
      for (auto &x : p) file << x << " ";
      file << std::endl;
    }

    file << c.Q.size() << std::endl;
    for (auto &p : c.Q) {
      for (auto &x : p) file << x << " ";
      file << std::endl;
    }
    file.close();
  }

  void writeToFile(NurbsSurface const &s,std::string const &file_name, int ulevel = 10, int vlevel=10)
  {
    using namespace vector_ops;

    // compute mesh
    std::vector<std::vector<int>> mesh;
    auto vid = [ulevel](int i, int j) { return i + j*(ulevel+1); };

    for (int j = 0; j < vlevel; j++)
    {
      for (int i = 0; i < ulevel; i++)
      {
        mesh.push_back({vid(i,j), vid(i+1,j), vid(i+1,j+1), vid(i,j+1)});
      }
    }

    // compute points on surface 
    std::vector<std::vector<double>> pts;
    auto du = (s.uknot.back() - s.uknot[0])/ulevel;
    auto dv = (s.vknot.back() - s.vknot[0])/vlevel;
    auto u = 0.0; auto v = 0.0;

    for (int j = 0; j <= vlevel; j++)
    {
      u = 0.0;
      for (int i = 0; i <= ulevel; i++)
      {
        pts.push_back(SurfacePoint(u,v,s));
        u += du;
      }
      v += dv;
    }

    // write file
    std::ofstream file;
    file.open(file_name);
    file << mesh.size() << std::endl; // number of points in file
    for (auto &el : mesh) {
      for (auto &v : el) file << v << " ";
      file << std::endl;
    }

    file << pts.size() << std::endl;
    for (auto &p : pts) {
      for (auto &x : p) file << x << " ";
      file << std::endl;
    }

    auto cpts = s.Q.size();
    file << cpts << std::endl;
    for (auto &p : s.Q)
    {
      for (auto &x : p) file << x << " ";
      file << std::endl;
    }
    file.close();
  }
}