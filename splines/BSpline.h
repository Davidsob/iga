#pragma once

#include "utils/VectorOperations.h"

#include "Algorithms.h"

#include <iostream>
#include <fstream>

#include <string>
#include <vector>


struct BSplineCurve
{
  using vector = std::vector<double>;
  using matrix = std::vector<std::vector<double>>; // need a point type

  size_t dim() const { return Q.empty() ? 0 : Q[0].size(); }
  int p; //polynomial order
  vector knot;
  matrix Q; 

  friend std::ostream operator<<(std::ostream const &os, BSplineCurve const &curve);
};

struct BSplineSurface
{
  using vector = std::vector<double>;
  using matrix = std::vector<std::vector<double>>; // need a point type

  size_t dim() const { return Q.empty() ? 0 : Q[0].size(); }
  int qid(int iu, int iv) const { return iu + iv*(uknot.size()-p-1); }
  int p, q; //polynomial order
  vector uknot, vknot; // knot vectors
  matrix Q; // cpts in vector form {c0j, c1j, cij...} 

  friend std::ostream operator<<(std::ostream const &os, BSplineSurface const &surf);
};

std::ostream & operator<<(std::ostream &os, BSplineCurve const &curve)
{
  using namespace vector_ops;
  os << "### BSpline-C ###" << std::endl;
  os << "p = " << curve.p << std::endl;
  os << "knot   = " << curve.knot << std::endl;
  os << "contol  = " << curve.Q << std::endl;
  return os;
}

std::ostream & operator<<(std::ostream &os, BSplineSurface const &surf)
{
  using namespace vector_ops;
  os << "### BSpline-S ###" << std::endl;
  os << "{p,q} = " << "{" << surf.p << "," << surf.q << "}" << std::endl;
  os << "uknot   = " << surf.uknot << std::endl;
  os << "vknot   = " << surf.vknot << std::endl;
  os << "contol  = " << surf.Q << std::endl;
  return os;
}

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

    point S(surf.dim(),0.0);
    std::vector<point> tmp(surf.q+1,S);

    for (int i = 0; i <= surf.q; i++)
    {
      for (int j = 0; j <= surf.p; j++)
      {
        auto a = uSpan-surf.p+j;
        auto b = vSpan-surf.q+i;
        auto idx = surf.qid(a,b);
        tmp[i] += Nu[j]*surf.Q[idx];
      }
    }

    for (int i = 0; i <= surf.q; i++)
    {
      S += Nv[i]*tmp[i];
    }
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

    point zero(surf.dim(),0);
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
          int idx = surf.qid(a,b);
          tmp[s] += Nu[k][r]*surf.Q[idx];
        }
        auto dd = std::min(order-k,dv);
        for (int l = 0; l <= dd; l++)
        {
          Skl[k][l] = zero;
          for (int s = 0; s <= surf.q; s++)
          {
            Skl[k][l] += Nv[l][s]*tmp[s];
          }
        }
      }
    }
    return Skl;
  }

  void writeToFile(BSplineCurve const &c,std::string const &file_name, int level = 20)
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

  void writeToFile(BSplineSurface const &s,std::string const &file_name, int ulevel = 10, int vlevel=10)
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

  void writeVectorData(std::vector<std::vector<double>> const &p,
                       std::vector<std::vector<double>> const &v,
                       std::string const &file_name, bool normalize=true, double scale = 1.0)
  {
    using namespace vector_ops;

    std::ofstream file;
    file.open(file_name);
    file << p.size() << std::endl;
    for (auto &x : p) {
      for (auto &xi : x) file << xi << " ";
      file << std::endl;
    }

    file << v.size() << std::endl;
    for (std::vector<double> const &x : v) {
      auto xn = !normalize ? x : scale*vector_ops::normalize(x);
      for (auto &xi : xn) file << xi << " ";
      file << std::endl;
    }
  }
}