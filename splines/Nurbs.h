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

  friend std::ostream operator<<(std::ostream const &os, NurbsCurve const &curve);
};

struct NurbsSurface 
{
  using vector = std::vector<double>;
  using matrix = std::vector<std::vector<double>>; // need a point type

  int qid(int iu, int iv) const { return iu + iv*(uknot.size()-p-1); }
  int p, q; //polynomial order
  vector uknot, vknot, weights; // knot vectors, weight vector
  matrix Q; // cpts in vector form {c0j, c1j, cij...} 

  friend std::ostream operator<<(std::ostream const &os, NurbsSurface const &surf);
};

struct NurbsSolid
{
  using vector = std::vector<double>;
  using matrix = std::vector<std::vector<double>>; // need a point type

  int qid(int iu, int iv) const { return iu + iv*(uknot.size()-p-1); }
  int p, q, r; //polynomial order
  vector uknot, vknot, wknot, weights; // knot vectors, weight vector
  matrix Q; // cpts in vector form {c0j, c1j, cij...} 

  friend std::ostream operator<<(std::ostream const &os, NurbsSolid const &solid);
};

std::ostream & operator<<(std::ostream &os, NurbsCurve const &curve)
{
  using namespace vector_ops;
  os << "### NURB-Curve ###" << std::endl;
  os << "p = " << curve.p << std::endl;
  os << "knot   = " << curve.knot << std::endl;
  os << "weights = " << curve.weights << std::endl;
  os << "contol  = " << curve.Q << std::endl;
  return os;
}

std::ostream & operator<<(std::ostream &os, NurbsSurface const &surf)
{
  using namespace vector_ops;
  os << "### NURB-Surface ###" << std::endl;
  os << "{p,q} = " << "{" << surf.p << "," << surf.q << "}" << std::endl;
  os << "uknot   = " << surf.uknot << std::endl;
  os << "vknot   = " << surf.vknot << std::endl;
  os << "weights = " << surf.weights << std::endl;
  os << "contol  = " << surf.Q << std::endl;
  return os;
}

std::ostream & operator<<(std::ostream &os, NurbsSolid const &solid)
{
  using namespace vector_ops;
  os << "### NURB-Solid ###" << std::endl;
  os << "{p,q,r} = " << "{" << solid.p << "," << solid.q << "," << solid.r << "}" << std::endl;
  os << "uknot   = " << solid.uknot << std::endl;
  os << "vknot   = " << solid.vknot << std::endl;
  os << "wknot   = " << solid.wknot << std::endl;
  os << "weights = " << solid.weights << std::endl;
  os << "contol  = " << solid.Q << std::endl;
  return os;
}

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
    b.q = s.q;
    b.uknot = s.uknot;
    b.vknot = s.vknot;
    weightedControlPoints(s,b.Q);
  }

  void weightedBSpline(NurbsSolid const &s, BSplineSolid &b)
  {
    b.p = s.p;
    b.q = s.q;
    b.r = s.r;
    b.uknot = s.uknot;
    b.vknot = s.vknot;
    b.wknot = s.wknot;
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
    auto Sw = spline_ops::SurfacePoint(u,v,b);
    Sw /= Sw.back(); Sw.pop_back();
    return Sw;
  }

  std::vector<double>
  SolidPoint(double u, double v, double w, NurbsSolid const &solid)
  {
    using namespace vector_ops;
    BSplineSolid b; weightedBSpline(solid,b);
    auto Sw = spline_ops::SolidPoint(u,v,w,b);
    Sw /= Sw.back(); Sw.pop_back();
    return Sw;
  }

  std::vector<double>
  SolidPoint2(double u, double v, double w, NurbsSolid const &solid)
  {
    using namespace vector_ops;
    BSplineSolid b; weightedBSpline(solid,b);
    auto Sw = spline_ops::SolidPoint2(u,v,w,b);
    Sw /= Sw.back(); Sw.pop_back();
    return Sw;
  }

  std::vector<std::vector<std::vector<double>>>
  SurfaceDerivatives(double u, double v, int order, NurbsSurface const &surf)
  {
    using namespace vector_ops;
    BSplineSurface b; weightedBSpline(surf,b);
    auto Aders = spline_ops::SurfaceDerivatives(u,v,order,b);

    size_t nk = Aders.size();
    size_t nl = Aders[0].size();
    using matrix = typename NurbsSurface::matrix;
    using point  = typename matrix::value_type;

    std::vector<matrix> Skl(nk, matrix(nl, point(3,0)));
     
    for (int k = 0; k <= order; k++)
    {
      for (int l = 0; l <= order-k; l++)
      {
        auto v = Aders[k][l];
        for (int j = 1; j <= l; j++)
        {
          v -= binomial(l,j)*Aders[0][j].back()*Skl[k][l-j];
        }

        for (int i = 1; i <= k; i++)
        {
          v -= binomial(k,i)*Aders[i][0].back()*Skl[k-i][l];
          decltype(v) v2(v.size(),0);
          for (int j = 1; j <= l; j++)
          {
            v2 += binomial(l,j)*Aders[i][j].back()*Skl[k-i][l-j];
          }
          v -= binomial(k,i)*v2;
        }
        v /= Aders[0][0].back(); v.pop_back();
        Skl[k][l] = v; 
      }
    }

    return Skl;
  }

  std::vector<double>
  SolidDerivative(double u, double v, double w, int order, int direction, NurbsSolid const &solid)
  {
    using namespace vector_ops;
    BSplineSolid b; weightedBSpline(solid,b);
    auto P = spline_ops::SolidPoint(u,v,w,b);
    auto Aders = spline_ops::SolidDerivative(u,v,w,order,direction,b);
    auto dP = (Aders - binomial(1,1)*Aders.back()*P)/P.back();
    dP.pop_back();
    return dP;
  }

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

  void writeToFile(NurbsSolid const &s,std::string const &file_name, int ulevel = 10, int vlevel=10, int wlevel=10)
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