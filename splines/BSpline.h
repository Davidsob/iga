#pragma once

#include "utils/VectorOperations.h"

#include "Algorithms.h"

#include <iostream>
#include <fstream>

#include <string>
#include <vector>

#include "GeometricObject.h"

struct BSplinePoint : public GeometricObject
{
  using vector = std::vector<double>;
  using matrix = std::vector<vector>; // need a point type

  explicit BSplinePoint() = default;
  explicit BSplinePoint(vector const x) : Q({x}) {}

  matrix const &coordinates() const override { return Q; }
  matrix Q; 

  friend std::ostream &operator<<(std::ostream &os, BSplinePoint const &point);
};

struct BSplineCurve : public GeometricObject
{
  using vector = std::vector<double>;
  using matrix = std::vector<vector>; // need a point type

  matrix const &coordinates() const override { return Q; }

  int p; //polynomial order
  vector knot;
  matrix Q; 

  friend std::ostream &operator<<(std::ostream &os, BSplineCurve const &curve);
};

struct BSplineSurface : public GeometricObject
{
  using vector = std::vector<double>;
  using matrix = std::vector<std::vector<double>>; // need a point type

  matrix const &coordinates() const override { return Q; }

  int qid(int iu, int iv) const { return iu + iv*(uknot.size()-p-1); }

  int p, q; //polynomial order
  vector uknot, vknot; // knot vectors
  matrix Q; // cpts in vector form {c0j, c1j, cij...} 

  friend std::ostream &operator<<(std::ostream &os, BSplineSurface const &surf);
};

struct BSplineSolid : public GeometricObject
{
  using vector = std::vector<double>;
  using matrix = std::vector<std::vector<double>>; // need a point type

  matrix const &coordinates() const override { return Q; }

  int qid(int a, int b, int c) const
  { 
    auto const cols = uknot.size()-p-1;
    auto const rows = vknot.size()-q-1;
    return a + b*(cols) + c*(cols*rows); 
  }

  int p, q, r; //polynomial order
  vector uknot, vknot, wknot; // knot vectors
  matrix Q; // cpts in vector form {c0j, c1j, cij...} 

  friend std::ostream &operator<<(std::ostream &os, BSplineSolid const &surf);
};

inline std::ostream & operator<<(std::ostream &os, BSplinePoint const &point)
{
  using namespace vector_ops;
  os << "### BSpline-Point ###" << std::endl;
  os << "contol  = " << point.Q << std::endl;
  return os;
}

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

    auto uspan = algo::FindSpan(u, surf.p, surf.uknot);
    auto vspan = algo::FindSpan(v, surf.q, surf.vknot);
    std::vector<double> Nu = algo::BasisFunctions(u,uspan,surf.p, surf.uknot);
    std::vector<double> Nv = algo::BasisFunctions(v,vspan,surf.q, surf.vknot);

    point S(surf.dim(),0.0);
    std::vector<point> tmp(surf.q+1,S);

    for (int i = 0; i <= surf.q; i++)
    {
      for (int j = 0; j <= surf.p; j++)
      {
        auto a = uspan-surf.p+j;
        auto b = vspan-surf.q+i;
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

    auto const uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    auto const vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    auto const wSpan = algo::FindSpan(w, solid.r, solid.wknot);
    std::vector<double> const Nu = algo::BasisFunctions(u,uSpan,solid.p, solid.uknot);
    std::vector<double> const Nv = algo::BasisFunctions(v,vSpan,solid.q, solid.vknot);
    std::vector<double> const Nw = algo::BasisFunctions(w,wSpan,solid.r, solid.wknot);

    point S(solid.dim(),0.0);

    for (int k = 0; k <= solid.r; k++)
    {
      auto c = k+wSpan-solid.r;
      for (int j = 0; j <= solid.q; j++)
      {
        auto b = j+vSpan-solid.q;
        for (int i = 0; i <= solid.p; i++)
        {
          auto a = i+uSpan-solid.p;
          auto idx = solid.qid(a,b,c);
          S += Nu[i]*Nv[j]*Nw[k]*solid.Q[idx];
        }
      }
    }

    return S;
  }

  std::vector<double>
  inline _SurfaceDerivativeWrtU(double u, double v, int order, BSplineSurface const &surf)
  {
    using matrix = typename BSplineSurface::matrix;
    using point  = typename matrix::value_type;

    auto const uSpan = algo::FindSpan(u, surf.p, surf.uknot);
    auto const vSpan = algo::FindSpan(v, surf.q, surf.vknot);

    auto const dNu = algo::BasisFunctionDerivatives(u,uSpan,surf.p,surf.uknot);
    auto const Nv = algo::BasisFunctions(v,vSpan,surf.q,surf.vknot);

    point dSdU(surf.dim(),0.0);

    for (int j = 0; j <= surf.q; j++)
    {
      auto b = j+vSpan-surf.q;
      for (int i = 0; i <= surf.p; i++)
      {
        auto a = i+uSpan-surf.p;
        auto idx = surf.qid(a,b);
        dSdU += dNu[order][i]*Nv[j]*surf.Q[idx];
      }
    }

    return dSdU;
  }

  std::vector<double>
  inline _SurfaceDerivativeWrtV(double u, double v, int order, BSplineSurface const &surf)
  {
    using matrix = typename BSplineSurface::matrix;
    using point  = typename matrix::value_type;

    auto const uSpan = algo::FindSpan(u, surf.p, surf.uknot);
    auto const vSpan = algo::FindSpan(v, surf.q, surf.vknot);

    auto const Nu = algo::BasisFunctions(u,uSpan,surf.p,surf.uknot);
    auto const dNv = algo::BasisFunctionDerivatives(v,vSpan,surf.q,surf.vknot);

    point dSdV(surf.dim(),0.0);

    for (int j = 0; j <= surf.q; j++)
    {
      auto b = j+vSpan-surf.q;
      for (int i = 0; i <= surf.p; i++)
      {
        auto a = i+uSpan-surf.p;
        auto idx = surf.qid(a,b);
        dSdV += Nu[i]*dNv[order][j]*surf.Q[idx];
      }
    }

    return dSdV;
  }

  std::vector<double>
  inline SurfaceDerivative(double u, double v, int order, int direction, BSplineSurface const &surf)
  {
    switch (direction)
    {
      case 1:  return _SurfaceDerivativeWrtV(u,v,order,surf); break;
      default: return _SurfaceDerivativeWrtU(u,v,order,surf); break;
    }
  }

  std::vector<double>
  inline _SolidDerivativeWrtU(double u, double v, double w, int order, BSplineSolid const &solid)
  {
    using matrix = typename BSplineSurface::matrix;
    using point  = typename matrix::value_type;

    auto const uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    auto const vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    auto const wSpan = algo::FindSpan(w, solid.r, solid.wknot);

    auto const dNu = algo::BasisFunctionDerivatives(u,uSpan,solid.p,solid.uknot);
    auto const Nv = algo::BasisFunctions(v,vSpan,solid.q,solid.vknot);
    auto const Nw = algo::BasisFunctions(w,wSpan,solid.r,solid.wknot);

    point dSdU(solid.dim(),0.0);

    for (int k = 0; k <= solid.r; k++)
    {
      auto c = k+wSpan-solid.r;
      for (int j = 0; j <= solid.q; j++)
      {
        auto b = j+vSpan-solid.q;
        for (int i = 0; i <= solid.p; i++)
        {
          auto a = i+uSpan-solid.p;
          auto idx = solid.qid(a,b,c);
          dSdU += dNu[order][i]*Nv[j]*Nw[k]*solid.Q[idx];
        }
      }
    }

    return dSdU;
  }

  std::vector<double>
  inline _SolidDerivativeWrtV(double u, double v, double w, int order, BSplineSolid const &solid)
  {
    using matrix = typename BSplineSurface::matrix;
    using point  = typename matrix::value_type;

    auto const uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    auto const vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    auto const wSpan = algo::FindSpan(w, solid.r, solid.wknot);

    auto const Nu = algo::BasisFunctions(u,uSpan,solid.p,solid.uknot);
    auto const dNv = algo::BasisFunctionDerivatives(v,vSpan,solid.q,solid.vknot);
    auto const Nw = algo::BasisFunctions(w,wSpan,solid.r,solid.wknot);

    point dSdV(solid.dim(),0.0);

    for (int k = 0; k <= solid.r; k++)
    {
      auto c = k+wSpan-solid.r;
      for (int j = 0; j <= solid.q; j++)
      {
        auto b = j+vSpan-solid.q;
        for (int i = 0; i <= solid.p; i++)
        {
          auto a = i+uSpan-solid.p;
          auto idx = solid.qid(a,b,c);
          dSdV += Nu[i]*dNv[order][j]*Nw[k]*solid.Q[idx];
        }
      }
    }

    return dSdV;
  }

  std::vector<double>
  inline _SolidDerivativeWrtW(double u, double v, double w, int order, BSplineSolid const &solid)
  {
    using matrix = typename BSplineSurface::matrix;
    using point  = typename matrix::value_type;

    auto const uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    auto const vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    auto const wSpan = algo::FindSpan(w, solid.r, solid.wknot);

    auto const Nu = algo::BasisFunctions(u,uSpan,solid.p,solid.uknot);
    auto const Nv = algo::BasisFunctions(v,vSpan,solid.q,solid.vknot);
    auto const dNw = algo::BasisFunctionDerivatives(w,wSpan,solid.r,solid.wknot);

    point dSdW(solid.dim(),0.0);

    for (int k = 0; k <= solid.r; k++)
    {
      auto c = k+wSpan-solid.r;
      for (int j = 0; j <= solid.q; j++)
      {
        auto b = j+vSpan-solid.q;
        for (int i = 0; i <= solid.p; i++)
        {
          auto a = i+uSpan-solid.p;
          auto idx = solid.qid(a,b,c);
          dSdW += Nu[i]*Nv[j]*dNw[order][k]*solid.Q[idx];
        }
      }
    }

    return dSdW;
  }

  std::vector<double>
  inline SolidDerivative(double u, double v, double w, int order, int direction, BSplineSolid const &solid)
  {
    switch (direction)
    {
      case 1:  return _SolidDerivativeWrtV(u,v,w,order,solid); break;
      case 2:  return _SolidDerivativeWrtW(u,v,w,order,solid); break;
      default: return _SolidDerivativeWrtU(u,v,w,order,solid); break;
    }
  }

  inline void _constUIsoCurve(double u, BSplineSurface const &surf, BSplineCurve &curve)
  {
    using namespace vector_ops;
    using point = typename BSplineSurface::matrix::value_type;

    auto uspan = algo::FindSpan(u, surf.p, surf.uknot);
    std::vector<double> Nu = algo::BasisFunctions(u,uspan,surf.p, surf.uknot);

    int const M = surf.vknot.size()-surf.q-1;
    for (int j = 0; j < M; j++)
    {
      point tmp(surf.dim(),0);
      for (int i = 0; i <= surf.p; i++)
      {
        auto a = uspan-surf.p+i;
        auto idx = surf.qid(a,j);
        tmp += Nu[i]*surf.Q[idx];
      }
      curve.Q.push_back(tmp);
    }

    curve.p = surf.q;
    curve.knot = surf.vknot; 
  }

  inline void _constVIsoCurve(double v, BSplineSurface const &surf, BSplineCurve &curve)
  {
    using namespace vector_ops;
    using point = typename BSplineSurface::matrix::value_type;

    auto vspan = algo::FindSpan(v, surf.q, surf.vknot);
    std::vector<double> Nv= algo::BasisFunctions(v,vspan,surf.q,surf.vknot);

    int const N = surf.uknot.size()-surf.p-1;
    for (int i = 0; i < N; i++)
    {
      point tmp(surf.dim(),0);
      for (int j = 0; j <= surf.q; j++)
      {
        auto b = vspan-surf.q+j;
        auto idx = surf.qid(i,b);
        tmp += Nv[j]*surf.Q[idx];
      }
      curve.Q.push_back(tmp);
    }

    curve.p = surf.p;
    curve.knot = surf.uknot; 
  }

  inline void getIsoCurve(double uv, size_t direction, BSplineSurface const &surf, BSplineCurve &curve)
  {
    switch(direction)
    {
      case 1 : _constVIsoCurve(uv,surf,curve); break;
      default: _constUIsoCurve(uv,surf,curve); break;
    }
  }

  inline void _constUIsoSurface(double u, BSplineSolid const &solid, BSplineSurface &surf)
  {
    using namespace vector_ops;

    using matrix = typename BSplineSurface::matrix;
    using point  = typename matrix::value_type;

    auto const uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    std::vector<double> const Nu = algo::BasisFunctions(u,uSpan,solid.p, solid.uknot);

    int const M = solid.vknot.size()-solid.q-1;
    int const O = solid.wknot.size()-solid.r-1;
    for (int j = 0; j < M; j++)
    {
      for (int k = 0; k < O; k++)
      {
        point tmp(solid.dim(),0);
        for (int i = 0; i <= solid.p; i++)
        {
          auto a = i+uSpan-solid.p;
          auto idx = solid.qid(a,j,k);
          tmp += Nu[i]*solid.Q[idx];
        }

        surf.Q.push_back(tmp);
      }
    }

    surf.p = solid.r;
    surf.q = solid.q;

    surf.uknot = solid.wknot;
    surf.vknot = solid.vknot;
  }

  inline void _constVIsoSurface(double v, BSplineSolid const &solid, BSplineSurface &surf)
  {
    using namespace vector_ops;

    using matrix = typename BSplineSurface::matrix;
    using point  = typename matrix::value_type;

    auto const vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    std::vector<double> const Nv = algo::BasisFunctions(v,vSpan,solid.q, solid.vknot);

    int const N = solid.uknot.size()-solid.p-1;
    int const O = solid.wknot.size()-solid.r-1;
    for (int k = 0; k < O; k++)
    {
      for (int i = 0; i < N; i++)
      {
        point tmp(solid.dim(),0);
        for (int j = 0; j <= solid.q; j++)
        {
          auto b = j+vSpan-solid.q;
          auto idx = solid.qid(i,b,k);
          tmp += Nv[j]*solid.Q[idx];
        }
        surf.Q.push_back(tmp);
      }
    }

    surf.p = solid.p;
    surf.q = solid.r;

    surf.uknot = solid.uknot;
    surf.vknot = solid.wknot;
  }

  inline void _constWIsoSurface(double w, BSplineSolid const &solid, BSplineSurface &surf)
  {
    using namespace vector_ops;

    using matrix = typename BSplineSurface::matrix;
    using point  = typename matrix::value_type;

    auto const wSpan = algo::FindSpan(w, solid.r, solid.wknot);
    std::vector<double> const Nw = algo::BasisFunctions(w,wSpan,solid.r, solid.wknot);

    int const N = solid.uknot.size()-solid.p-1;
    int const M = solid.vknot.size()-solid.q-1;
    for (int j = 0; j < M; j++)
    {
      for (int i = 0; i < N; i++)
      {
        point tmp(solid.dim(),0);
        for (int k = 0; k <= solid.r; k++)
        {
          auto c = k+wSpan-solid.r;
          auto idx = solid.qid(i,j,c);
          tmp += Nw[k]*solid.Q[idx];
        }
        surf.Q.push_back(tmp);
      }
    }

    surf.p = solid.p;
    surf.q = solid.q;

    surf.uknot = solid.uknot;
    surf.vknot = solid.vknot;
  }

  inline void getIsoSurface(double uvw, size_t direction, BSplineSolid const &solid, BSplineSurface &surface)
  {
    switch(direction)
    {
      case 1 : _constVIsoSurface (uvw,solid,surface); break;
      case 2 : _constWIsoSurface (uvw,solid,surface); break;
      default: _constUIsoSurface(uvw,solid,surface); break;
    }
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
    file.close();
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