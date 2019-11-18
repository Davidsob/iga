#include "SplineModifiers.h"

#include "Algorithms.h"
#include "BSpline.h"
#include "Nurbs.h"
#include "utils/VectorOperations.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace spline_ops
{
  void _insertKnot(double u, BSplineCurve &shape)
  {
    auto const &p = shape.p;
    auto const &U = shape.knot;
    auto const &Pw = shape.Q;
    auto const  np = Pw.size();
    auto const   k = algo::FindSpan(u,p,U);

    std::vector<std::vector<double>> Qw;
    for (size_t i = 0; i <= k-p; i++) Qw.push_back(Pw[i]);
    for (size_t i = k-p+1; i <= k; i++)
    {
      auto alf = (u-U[i])/(U[i+p]-U[i]);
      Qw.push_back(alf*Pw[i] + (1.0-alf)*Pw[i-1]);
    }
    for (size_t i = k+1; i <= np; i++) Qw.push_back(Pw[i-1]);

    shape.knot.insert(shape.knot.begin()+k+1,u);
    shape.Q = Qw;
  }

  void insertKnot(double u, int r, BSplineCurve &shape)
  {
    if (r == 0) return;
    if (shape.knot.empty() || shape.Q.empty()) return;

    int const  k = algo::FindSpan(u,shape.p,shape.knot);
    int const  s = algo::SpanMultiplicity(k,shape.knot);

    if ((s+r) > shape.p) r = shape.p-s; // adjust number of repititions if required 

    while (r-- > 0)
    {
      _insertKnot(u,shape);
    }
  }

  void _insertKnotU(double u, BSplineSurface &shape)
  {
    auto const &p = shape.p;
    auto const &q = shape.q;
    auto const &U = shape.uknot;
    auto const &V = shape.vknot;
    auto const &Pw = shape.Q;
    auto const  k  = algo::FindSpan(u,p,U);

    // figure out how many rows and columns of points we are dealing with
    auto const nu = U.size()-p-1;
    auto const nv = V.size()-q-1;
    auto idx = [&nu](auto const iu, auto const iv)
    {
      return iu + (nu+1)*iv; 
    };

    std::vector<double> zero(shape.dim(),0);
    std::vector<std::vector<double>> Qw((nu+1)*nv,zero);

    for (size_t j = 0; j < nv; j++)
    {
      for (size_t i = 0; i <= k-p; i++) Qw[idx(i,j)] = Pw[shape.qid(i,j)];
      for (size_t i = k-p+1; i <= k; i++)
      {
        auto a = shape.qid(i,j);    // ith point
        auto b = shape.qid(i-1,j);  // i-1th point

        auto alf = (u-U[i])/(U[i+p]-U[i]);
        Qw[idx(i,j)] = (alf*Pw[a] + (1.0-alf)*Pw[b]);
      }
      for (size_t i = k+1; i <= nu; i++) Qw[idx(i,j)] = Pw[shape.qid(i-1,j)];
    }

    shape.uknot.insert(shape.uknot.begin()+k+1,u);
    shape.Q = Qw;
  }

  void _insertKnotV(double v, BSplineSurface &shape)
  {
    auto const &p = shape.p;
    auto const &q = shape.q;
    auto const &U = shape.uknot;
    auto const &V = shape.vknot;
    auto const &Pw = shape.Q;
    auto const  k  = algo::FindSpan(v,q,V);

    // figure out how many rows and columns of points we are dealing with
    auto const nu = U.size()-p-1;
    auto const nv = V.size()-q-1;
    auto idx = [&nu](auto const iu, auto const iv)
    {
      return iu + (nu)*iv; 
    };

    std::vector<double> zero(shape.dim(),0);
    std::vector<std::vector<double>> Qw((nu)*(nv+1),zero);

    for (size_t i = 0; i < nu; i++)
    {
      for (size_t j = 0; j <= k-q; j++) Qw[idx(i,j)] = Pw[shape.qid(i,j)];
      for (size_t j = k-q+1; j <= k; j++)
      {
        auto a = shape.qid(i,j);    // ith point
        auto b = shape.qid(i,j-1);  // i-1th point

        auto alf = (v-V[j])/(V[j+q]-V[j]);
        Qw[idx(i,j)] = (alf*Pw[a] + (1.0-alf)*Pw[b]);
      }
      for (size_t j = k+1; j <= nv; j++) Qw[idx(i,j)] = Pw[shape.qid(i,j-1)];
    }

    shape.vknot.insert(shape.vknot.begin()+k+1,v);
    shape.Q = Qw;
  }

  void insertKnot(double uv, int r, int direction, BSplineSurface &shape)
  {
    if (r == 0) return;
    if (shape.Q.empty()) return;
    if (direction < 0 || direction > 1) return;

    if (direction == 0)
    {
      if (shape.uknot.empty()) return;

      int const  k = algo::FindSpan(uv,shape.p,shape.uknot);
      int const  s = algo::SpanMultiplicity(k,shape.uknot);

      if ((s+r) > shape.p) r = shape.p-s; // adjust number of repititions if required 

      while (r-- > 0)
      {
        _insertKnotU(uv,shape);
      }
    }
    else
    {
      if (shape.vknot.empty()) return;

      int const  k = algo::FindSpan(uv,shape.q,shape.vknot);
      int const  s = algo::SpanMultiplicity(k,shape.vknot);

      if ((s+r) > shape.q) r = shape.q-s; // adjust number of repititions if required 

      while (r-- > 0)
      {
        _insertKnotV(uv,shape);
      }
    }
  }

  void _insertKnotU(double u, BSplineSolid &shape)
  {
    auto const &p = shape.p;
    auto const &q = shape.q;
    auto const &r = shape.r;
    auto const &U = shape.uknot;
    auto const &V = shape.vknot;
    auto const &W = shape.wknot;
    auto const &Pw = shape.Q;
    auto const  s  = algo::FindSpan(u,p,U);

    // figure out how many rows and columns of points we are dealing with
    auto const nu = U.size()-p-1;
    auto const nv = V.size()-q-1;
    auto const nw = W.size()-r-1;
    auto idx = [&nu, &nv](auto const iu, auto const iv, auto const iw)
    {
      return iu + (nu+1)*iv + iw*(nu+1)*nv; 
    };

    std::vector<double> zero(shape.dim(),0);
    std::vector<std::vector<double>> Qw((nu+1)*nv*nw,zero);

    for (size_t k = 0; k < nw; k++)
    {
      for (size_t j = 0; j < nv; j++)
      {
        for (size_t i = 0; i <= s-p; i++) Qw[idx(i,j,k)] = Pw[shape.qid(i,j,k)];
        for (size_t i = s-p+1; i <= s; i++)
        {
          auto a = shape.qid(i,j,k);    // ith point
          auto b = shape.qid(i-1,j,k);  // i-1th point

          auto alf = (u-U[i])/(U[i+p]-U[i]);
          Qw[idx(i,j,k)] = (alf*Pw[a] + (1.0-alf)*Pw[b]);
        }
        for (size_t i = s+1; i <= nu; i++) Qw[idx(i,j,k)] = Pw[shape.qid(i-1,j,k)];
      }
    } 

    shape.uknot.insert(shape.uknot.begin()+s+1,u);
    shape.Q = Qw;
  }

  void _insertKnotV(double v, BSplineSolid &shape)
  {
    auto const &p = shape.p;
    auto const &q = shape.q;
    auto const &r = shape.r;
    auto const &U = shape.uknot;
    auto const &V = shape.vknot;
    auto const &W = shape.wknot;
    auto const &Pw = shape.Q;
    auto const  s  = algo::FindSpan(v,q,V);

    // figure out how many rows and columns of points we are dealing with
    auto const nu = U.size()-p-1;
    auto const nv = V.size()-q-1;
    auto const nw = W.size()-r-1;
    auto idx = [&nu, &nv](auto const iu, auto const iv, auto const iw)
    {
      return iu + (nu)*iv + iw*(nu)*(nv+1); 
    };

    std::vector<double> zero(shape.dim(),0);
    std::vector<std::vector<double>> Qw(nu*(nv+1)*nw,zero);

    for (size_t i = 0; i < nu; i++)
    {
      for (size_t k = 0; k < nw; k++)
      {
        for (size_t j = 0; j <= s-q; j++) Qw[idx(i,j,k)] = Pw[shape.qid(i,j,k)];
        for (size_t j = s-q+1; j <= s; j++)
        {
          auto a = shape.qid(i,j,k);    // jth pojnt
          auto b = shape.qid(i,j-1,k);  // j-1th pojnt

          auto alf = (v-V[j])/(V[j+q]-V[j]);
          Qw[idx(i,j,k)] = (alf*Pw[a] + (1.0-alf)*Pw[b]);
        }
        for (size_t j = s+1; j <= nv; j++) Qw[idx(i,j,k)] = Pw[shape.qid(i,j-1,k)];
      }
    } 

    shape.vknot.insert(shape.vknot.begin()+s+1,v);
    shape.Q = Qw;
  }

  void _insertKnotW(double w, BSplineSolid &shape)
  {
    auto const &p = shape.p;
    auto const &q = shape.q;
    auto const &r = shape.r;
    auto const &U = shape.uknot;
    auto const &V = shape.vknot;
    auto const &W = shape.wknot;
    auto const &Pw = shape.Q;
    auto const  s  = algo::FindSpan(w,r,W);

    // figure out how many rows and columns of points we are dealing with
    auto const nu = U.size()-p-1;
    auto const nv = V.size()-q-1;
    auto const nw = W.size()-r-1;
    auto idx = [&nu, &nv](auto const iu, auto const iv, auto const iw)
    {
      return iu + nu*iv + iw*nu*nv; 
    };

    std::vector<double> zero(shape.dim(),0);
    std::vector<std::vector<double>> Qw(nu*nv*(nw+1),zero);

    for (size_t j = 0; j < nv; j++)
    {
      for (size_t i = 0; i < nu; i++)
      {
        for (size_t k = 0; k <= s-r; k++) Qw[idx(i,j,k)] = Pw[shape.qid(i,j,k)];
        for (size_t k = s-r+1; k <= s; k++)
        {
          auto a = shape.qid(i,j,k);    // kth poknt
          auto b = shape.qid(i,j,k-1);  // k-1th poknt

          auto alf = (w-W[k])/(W[k+r]-W[k]);
          Qw[idx(i,j,k)] = (alf*Pw[a] + (1.0-alf)*Pw[b]);
        }
        for (size_t k = s+1; k <= nw; k++) Qw[idx(i,j,k)] = Pw[shape.qid(i,j,k-1)];
      }
    } 

    shape.wknot.insert(shape.wknot.begin()+s+1,w);
    shape.Q = Qw;
  }

  void insertKnot(double uvw, int r, int direction, BSplineSolid &shape)
  {
    if (r == 0) return;
    if (shape.Q.empty()) return;
    if (direction < 0 || direction > 2) return;

    if (direction == 0)
    {
      if (shape.uknot.empty()) return;

      int const  k = algo::FindSpan(uvw,shape.p,shape.uknot);
      int const  s = algo::SpanMultiplicity(k,shape.uknot);

      if ((s+r) > shape.p) r = shape.p-s; // adjust number of repititions if required 

      while (r-- > 0)
      {
        _insertKnotU(uvw,shape);
      }
    }
    else if (direction == 1)
    {
      if (shape.vknot.empty()) return;

      int const  k = algo::FindSpan(uvw,shape.q,shape.vknot);
      int const  s = algo::SpanMultiplicity(k,shape.vknot);

      if ((s+r) > shape.q) r = shape.q-s; // adjust number of repititions if required 

      while (r-- > 0)
      {
        _insertKnotV(uvw,shape);
      }
    } 
    else
    {
      if (shape.wknot.empty()) return;

      int const  k = algo::FindSpan(uvw,shape.r,shape.wknot);
      int const  s = algo::SpanMultiplicity(k,shape.wknot);

      if ((s+r) > shape.r) r = shape.r-s; // adjust number of repititions if required 

      while (r-- > 0)
      {
        _insertKnotW(uvw,shape);
      }
    }
  }

  void insertKnot(double u, size_t rep, NurbsCurve &curve)
  {
    BSplineCurve b; spline_ops::weightedBSpline(curve,b);
    insertKnot(u,rep,b);
    curve.knot = b.knot;
    spline_ops::extractWeights(b.Q,curve.Q,curve.weights);
  }

  void insertKnot(double uv, size_t rep, size_t direction, NurbsSurface &surface)
  {
    BSplineSurface b; spline_ops::weightedBSpline(surface,b);
    insertKnot(uv,rep,direction,b);
    surface.uknot = b.uknot;
    surface.vknot = b.vknot;
    spline_ops::extractWeights(b.Q,surface.Q,surface.weights);
  }

  void insertKnot(double uvw, size_t rep, size_t direction, NurbsSolid &solid)
  {
    BSplineSolid b; spline_ops::weightedBSpline(solid,b);
    insertKnot(uvw,rep,direction,b);
    solid.uknot = b.uknot;
    solid.vknot = b.vknot;
    solid.wknot = b.wknot;
    spline_ops::extractWeights(b.Q,solid.Q,solid.weights);
  }

  std::vector<double> knotMidpoints(std::vector<double> const &knot)
  {
    std::vector<double> midpoints;

    for (size_t i = 1; i < knot.size(); i++)
    {
      if (!algo::equal(knot[i],knot[i-1]))
        midpoints.push_back(0.5*(knot[i]+knot[i-1]));
    }

    return midpoints;
  }


  template<typename Curve>
  void _curveRefinement(size_t level, Curve& curve)
  {
    while (level-- > 0)
    {
      for (auto &mp : knotMidpoints(curve.knot))
      {
        insertKnot(mp,1,curve);
      }
    }
  }

  template<typename Shape>
  void _anisotropicRefinement(size_t level, size_t direction, Shape& shape, std::vector<double> &knot)
  {
    // the selector provides the correct knot vector for the specified direction
    while (level-- > 0)
    {
      for (auto &mp : knotMidpoints(knot))
      {
        insertKnot(mp,1,direction,shape);
      }
    }
  }

  template<>
  void midpointRefinement(size_t level, BSplineCurve& curve)
  {
    _curveRefinement(level,curve);
  }

  template<>
  void midpointRefinement(size_t level, NurbsCurve& curve)
  {
    _curveRefinement(level,curve);
  }

  template<>
  void midpointRefinement(size_t level, size_t direction, BSplineSurface &shape)
  {
    if (direction == 0)
    {
      _anisotropicRefinement(level,direction,shape,shape.uknot);
    } else {
      _anisotropicRefinement(level,direction,shape,shape.vknot);
    }
  }

  template<>
  void midpointRefinement(size_t level, size_t direction, NurbsSurface &shape)
  {
    if (direction == 0)
    {
      _anisotropicRefinement(level,direction,shape,shape.uknot);
    } else {
      _anisotropicRefinement(level,direction,shape,shape.vknot);
    }
  }

  template<>
  void midpointRefinement(size_t level, size_t direction, BSplineSolid &shape)
  {
    if (direction == 0)
    {
      _anisotropicRefinement(level,direction,shape,shape.uknot);
    } else {
      _anisotropicRefinement(level,direction,shape,shape.vknot);
    }
  }

  template<>
  void midpointRefinement(size_t level, size_t direction, NurbsSolid &shape)
  {
    if (direction == 0)
    {
      _anisotropicRefinement(level,direction,shape,shape.uknot);
    } else if (direction == 1) {
      _anisotropicRefinement(level,direction,shape,shape.vknot);
    } else {
      _anisotropicRefinement(level,direction,shape,shape.wknot);
    }
  }
}
