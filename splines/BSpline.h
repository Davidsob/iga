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

struct BSplineSolid
{
  using vector = std::vector<double>;
  using matrix = std::vector<std::vector<double>>; // need a point type

  size_t dim() const { return Q.empty() ? 0 : Q[0].size(); }
  int qid(int a, int b, int c) const
  { 
    auto cols = uknot.size()-p-1;
    auto rows = vknot.size()-q-1;
    return a + b*(cols) + c*(cols*rows); 
  }
  int p, q, r; //polynomial order
  vector uknot, vknot, wknot; // knot vectors
  matrix Q; // cpts in vector form {c0j, c1j, cij...} 

  friend std::ostream operator<<(std::ostream const &os, BSplineSolid const &surf);
};

inline std::ostream & operator<<(std::ostream &os, BSplineCurve const &curve)
{
  using namespace vector_ops;
  os << "### BSpline-Curve ###" << std::endl;
  os << "p = " << curve.p << std::endl;
  os << "knot   = " << curve.knot << std::endl;
  os << "contol  = " << curve.Q << std::endl;
  return os;
}

inline std::ostream & operator<<(std::ostream &os, BSplineSurface const &surf)
{
  using namespace vector_ops;
  os << "### BSpline-Surface ###" << std::endl;
  os << "{p,q} = " << "{" << surf.p << "," << surf.q << "}" << std::endl;
  os << "uknot   = " << surf.uknot << std::endl;
  os << "vknot   = " << surf.vknot << std::endl;
  os << "contol  = " << surf.Q << std::endl;
  return os;
}

inline std::ostream & operator<<(std::ostream &os, BSplineSolid const &solid)
{
  using namespace vector_ops;
  os << "### BSpline-Solid ###" << std::endl;
  os << "{p,q,r} = " << "{" << solid.p << "," << solid.q << "," << solid.r << "}" << std::endl;
  os << "uknot   = " << solid.uknot << std::endl;
  os << "vknot   = " << solid.vknot << std::endl;
  os << "wknot   = " << solid.wknot << std::endl;
  os << "contol  = " << solid.Q << std::endl;
  return os;
}

namespace spline_ops
{
  std::vector<double>
  inline CurvePoint(double u, BSplineCurve const &curve)
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
  inline CurveDerivatives(double u, int order, BSplineCurve const &curve)
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
  inline SurfacePoint(double u, double v, BSplineSurface const &surf)
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
  inline SurfaceDerivatives(double u, double v, int order, BSplineSurface const &surf)
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

  std::vector<double>
  inline SolidPoint(double u, double v, double w, BSplineSolid const &solid)
  {
    using namespace vector_ops;

    using matrix = typename BSplineSurface::matrix;
    using point  = typename matrix::value_type;

    auto uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    auto vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    auto wSpan = algo::FindSpan(w, solid.r, solid.wknot);
    std::vector<double> Nu = algo::BasisFunctions(u,uSpan,solid.p, solid.uknot);
    std::vector<double> Nv = algo::BasisFunctions(v,vSpan,solid.q, solid.vknot);
    std::vector<double> Nw = algo::BasisFunctions(w,wSpan,solid.r, solid.wknot);

    point S(solid.dim(),0.0);
    std::vector<matrix> tmp(solid.q+1,matrix(solid.r+1,S));

    for (int k = 0; k <= solid.r; k++)
    {
      auto c = wSpan-solid.r+k;
      for (int j = 0; j <= solid.q; j++)
      {
        auto b = vSpan-solid.q+j;
        for (int i = 0; i <= solid.p; i++)
        {
          auto a = uSpan-solid.p+i;
          auto idx = solid.qid(a,b,c);
          tmp[j][k] += Nu[i]*solid.Q[idx];
        }
      }
    }

    for (int j = 0; j <= solid.r; j++)
    {
      for (int i = 0; i <= solid.q; i++)
      {
        S += Nw[j]*Nv[i]*tmp[i][j];
      }
    }
    return S;
  }

  std::vector<double>
  inline SolidPoint2(double u, double v, double w, BSplineSolid const &solid)
  {
    using namespace vector_ops;

    using matrix = typename BSplineSurface::matrix;
    using point  = typename matrix::value_type;

    point S(solid.dim(),0.0);

    int n = solid.uknot.size()-solid.p-2;
    int m = solid.vknot.size()-solid.q-2;
    int o = solid.wknot.size()-solid.r-2;

    for (int k = 0; k <= o; k++)
    {
      auto Nw = algo::BasisFunction(w,k,solid.r,solid.wknot);
      for (int j = 0; j <= m; j++)
      {
        auto Nv = algo::BasisFunction(v,j,solid.q,solid.vknot);
        for (int i = 0; i <= n; i++)
        {
          auto idx = solid.qid(i,j,k);
          auto Nu = algo::BasisFunction(u,i,solid.p,solid.uknot);
          S += Nu*Nv*Nw*solid.Q[idx];
        }
      }
    }

    return S;
  }

  std::vector<double>
  inline SolidDerivative(double u, double v, double w, int order, int direction, BSplineSolid const &solid)
  {
    using namespace vector_ops;

    using matrix = typename BSplineSurface::matrix;
    using point  = typename matrix::value_type;

    auto basis = [](double u, double i, double p, auto const &knot)
    {
      return algo::BasisFunction(u,i,p,knot);
    };

    auto derivative = [order](double u, double i, double p, auto const &knot)
    {
      return algo::BasisFunctionDerivative(u,i,p,knot)[order];
    };

    point dS(solid.dim(),0.0);
    int n = solid.uknot.size()-solid.p-2;
    int m = solid.vknot.size()-solid.q-2;
    int o = solid.wknot.size()-solid.r-2;

    switch(direction)
    {
      case 1:
      {
        for (int k = 0; k <= o; k++)
        {
          auto Nw = basis(w,k,solid.r,solid.wknot);
          for (int j = 0; j <= m; j++)
          {
            auto dNv = derivative(v,j,solid.q,solid.vknot);
            for (int i = 0; i <= n; i++)
            {
              auto idx = solid.qid(i,j,k);
              auto Nu = basis(u,i,solid.p,solid.uknot);
              dS += Nu*dNv*Nw*solid.Q[idx];
            }
          }
        }
        break;
      }

      case 2:
      {
        for (int k = 0; k <= o; k++)
        {
          auto dNw = derivative(w,k,solid.r,solid.wknot);
          for (int j = 0; j <= m; j++)
          {
            auto Nv = basis(v,j,solid.q,solid.vknot);
            for (int i = 0; i <= n; i++)
            {
              auto idx = solid.qid(i,j,k);
              auto Nu = basis(u,i,solid.p,solid.uknot);
              dS += Nu*Nv*dNw*solid.Q[idx];
            }
          }
        }
        break;
      }

     default:
      {
        for (int k = 0; k <= o; k++)
        {
          auto Nw = basis(w,k,solid.r,solid.wknot);
          for (int j = 0; j <= m; j++)
          {
            auto Nv = basis(v,j,solid.q,solid.vknot);
            for (int i = 0; i <= n; i++)
            {
              auto idx = solid.qid(i,j,k);
              auto dNu = derivative(u,i,solid.p,solid.uknot);
              dS += dNu*Nv*Nw*solid.Q[idx];
            }
          }
        }
        break;
      }
    }

    return dS;
  }

  inline void writeToFile(BSplineCurve const &c,std::string const &file_name, int level = 20)
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

  inline void writeToFile(BSplineSurface const &s,std::string const &file_name, int ulevel = 10, int vlevel=10)
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

  inline void writeVectorData(std::vector<std::vector<double>> const &p,
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

  inline void writeToFile(BSplineSolid const &s,std::string const &file_name, int ulevel = 10, int vlevel=10, int wlevel=10)
  {
    using namespace vector_ops;

    // compute mesh
    std::vector<std::vector<int>> mesh;
    auto cols = ulevel+1;
    auto rows = vlevel+1;
    auto vid = [cols,rows](int i, int j, int k) { return i + j*(cols) + k*(cols*rows); };

    // plot surfaces!
    // x-
    for (int j = 0; j < wlevel; j++)
    {
      for (int i = 0; i < vlevel; i++)
      {
        mesh.push_back({vid(0,i,j), vid(0,i+1,j), vid(0,i+1,j+1), vid(0,i,j+1)});
      }
    }
    // x+
    for (int j = 0; j < wlevel; j++)
    {
      for (int i = 0; i < vlevel; i++)
      {
        mesh.push_back({vid(ulevel,i,j), vid(ulevel,i+1,j), vid(ulevel,i+1,j+1), vid(ulevel,i,j+1)});
      }
    }
    // y-
    for (int j = 0; j < ulevel; j++)
    {
      for (int i = 0; i < wlevel; i++)
      {
        mesh.push_back({vid(j,0,i), vid(j,0,i+1), vid(j+1,0,i+1), vid(j+1,0,i)});
      }
    }
    // y+
    for (int j = 0; j < ulevel; j++)
    {
      for (int i = 0; i < wlevel; i++)
      {
        mesh.push_back({vid(j,vlevel,i), vid(j,vlevel,i+1), vid(j+1,vlevel,i+1), vid(j+1,vlevel,i)});
      }
    }
    // z-
    for (int j = 0; j < vlevel; j++)
    {
      for (int i = 0; i < ulevel; i++)
      {
        mesh.push_back({vid(i,j,0), vid(i+1,j,0), vid(i+1,j+1,0), vid(i,j+1,0)});
      }
    }
    // z+
    for (int j = 0; j < vlevel; j++)
    {
      for (int i = 0; i < ulevel; i++)
      {
        mesh.push_back({vid(i,j,wlevel), vid(i+1,j,wlevel), vid(i+1,j+1,wlevel), vid(i,j+1,wlevel)});
      }
    }

    // compute points on surface 
    std::vector<std::vector<double>> pts;
    auto du = (s.uknot.back() - s.uknot[0])/ulevel;
    auto dv = (s.vknot.back() - s.vknot[0])/vlevel;
    auto dw = (s.wknot.back() - s.wknot[0])/wlevel;
    auto u = 0.0; auto v = 0.0; auto w = 0.0;

    for (int k = 0; k <= wlevel; k++)
    {
      v = 0.0;
      for (int j = 0; j <= vlevel; j++)
      {
        u = 0.0;
        for (int i = 0; i <= ulevel; i++)
        {
          pts.push_back(SolidPoint(u,v,w,s));
          u += du;
        }
        v += dv;
      }
      w += dw;
    }

    // write file
    std::ofstream file;
    file.open(file_name);
    file << mesh.size() << std::endl; // number of surface elements
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