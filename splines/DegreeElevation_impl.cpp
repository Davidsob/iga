// #include "SplineModifiers.h"

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
  void elevate(size_t rep, BSplineCurve &curve)
  {
    if (rep == 0) return;

    auto const p = curve.p;
    auto const &Pw = curve.Q;
    auto const &U = curve.knot;
    auto const n = U.size()-p-1;

    int const m = n+p;
    int const t = rep;
    int const ph = p+t;

    // compute bezier degree elevation coefficeints 
    std::vector<double> zero(curve.dim(),0);
    // std::vector<double> Uh(U.size()+2*t,0.0);
    std::vector<double> Uh;
    std::vector<std::vector<double>> Qw;
    // std::vector<std::vector<double>> Qw(m+2,zero);
    std::vector<std::vector<double>> bezalfs(p+t+1, std::vector<double>(p+1,0));
    std::vector<std::vector<double>> bpts(p+1, zero);
    std::vector<std::vector<double>> ebpts(p+t+1 , zero);
    std::vector<std::vector<double>> nextBpts(p-1, zero);
    std::vector<double> alphas(p-1,0);


    // compute bez elevation coefss
    int const ph2 = ph/2;
    bezalfs[0][0] = bezalfs[ph][p] = 1.0;
    for (int i = 1; i <= ph2 ; i++)
    {
      auto inv = 1.0/algo::binomial(ph,i);
      auto mpi = std::min(p,i);
      for (int j = std::max(0,i-t); j <= mpi; j++)
      {
        bezalfs[i][j] = inv*algo::binomial(p,j)*algo::binomial(t,i-j);
      }
    }

    for (int i = ph2+1; i <= ph-1 ; i++)
    {
      auto mpi = std::min(p,i);
      for (int j = std::max(0,i-t); j <= mpi; j++)
      {
        bezalfs[i][j] = bezalfs[ph-i][p-j];
      }
    }

    auto mh = ph;
    auto ua = U[0];
    for (int i = 0; i <= ph; i++) Uh.push_back(ua);
    auto kind = ph+1;
    Qw.push_back(Pw[0]);
    auto cind = 1;

    // initialize first bez segment 
    int lbz = 1; int r = -1;
    for (int i = 0; i <= p; i++) bpts[i] = Pw[i];
    auto a = p; auto b = p+1;
    while(b < m) // big loop through knot
    {
      auto i = b;
      while(b < m && algo::equal(U[b],U[b+1]))
      {
        b = b+1;
      }
      auto mul = b-i+1;
      mh = mh + mul+t;
      auto ub = U[b];
      auto oldr = r;
      r = p-mul; // insert knot b r times
      lbz = (oldr > 0) ? (oldr+2)/2 : 1;
      auto rbz = (r > 0) ? ph-(r+1)/2 : ph;

      if (r > 0) // begin insert knot
      {
        auto numer = ub-ua;
        for (int k = p; k > mul; k--)
          alphas[k-mul-1] = numer/(U[a+k]-ua);
        for (int j = 1; j <= r; j++)
        {
          auto save = r-j;
          auto s = mul+j;
          for (int k = p; k >= s; k--)
          {
            auto alf = alphas[k-s];
            bpts[k] = alf*bpts[k] + (1.0-alf)*bpts[k-1];
          }
          nextBpts[save] = bpts[p];
        }
      } // end insert knot

      auto RBZ = (r > 0) ? rbz+1 : rbz;  // degree elevate bez curve
      for (i = lbz; i <= RBZ; i++)
      {
        ebpts[i] = zero;
        auto mpi = std::min(p,i);
        for (int j = std::max(0,i-t); j <= mpi; j++)
        {
          ebpts[i] += bezalfs[i][j]*bpts[j];
        } // end of degree elevation of curv
      }

      if (oldr > 1) // begin removing knot u = U[a]
      {
        // must remove knot u= U[a] oldr times
        auto alfj = 1.0;
        if (oldr > 2)
        {
          alfj = (ua - Uh[kind-1])/(ub-Uh[kind-1]);
        }

        auto first = kind-2;
        auto last = kind;
        for (int tr = 1; tr < oldr; tr++)
        {
          i = first;
          auto j = last;
          auto kj = j-kind+1;
          while (j-1 > tr) // compute new control points for one removal step
          {
            if (i < cind)
            {
              auto alf = (ua-Uh[i])/(ub-Uh[i]);
              Qw[i] = (Qw[i] - (1.0-alf)*Qw[i-1])/alf;
            }

            if (kj >= lbz)
            {
              ebpts[kj] -= alfj*ebpts[kj+1]/(1.0-alfj);
            }

            i++; j--; kj--;
          }
          first--; last++;
        }
      } // end of removoing knot

      if (a != p) // load knot ua
      {
        for (i = 0; i < ph-oldr; i++)
        {
          // Uh[kind] = ua;
          Uh.push_back(ua);
          kind++;
        }
      }

      for (int j = lbz; j <= rbz; j++)
      {
        Qw.push_back(ebpts[j]);
        cind++;
      }

      if (b < m) // set up for next pass through loop
      {
        lbz = 1;
        for (int j = 0; j < r; j++)  bpts[j] = nextBpts[j];
        for (int j = r; j <= p; j++) bpts[j] = Pw[b-p+j];
        a = b;
        b++;
        ua = ub;
      } else // end knot
      {
        for (i = 0; i <= ph; i++) {
          Uh.push_back(ub);
        }
      } 
    } // end of while(b<m)

    // set new curve
    curve.p += t;
    curve.knot = Uh;
    curve.Q = Qw;
  }

  void _elevateU(int t, BSplineSurface &shape)
  {
    // elevate iso-curves of constant v and then
    // reconstruct the elevated surface patch using
    // the elevated curves -- easy peasy.
    auto const n = shape.uknot.size()-shape.p-1;
    auto const m = shape.vknot.size()-shape.q-1;
    // load v-const iso-curves and then elevate
    size_t j = 0;
    std::vector<BSplineCurve> curves(m,BSplineCurve());
    for (BSplineCurve &curve : curves)
    {
      curve.p = shape.p;
      curve.knot = shape.uknot;
      for (size_t i = 0; i < n; i++) curve.Q.push_back(shape.Q[shape.qid(i,j)]);
      elevate(t,curve);
      j++;
    }
    // assemble the elevated surface
    shape.p += t;
    shape.uknot = curves[0].knot;
    shape.Q.clear();

    for (BSplineCurve const &curve : curves)
    {
      shape.Q.insert(shape.Q.end(),curve.Q.begin(),curve.Q.end());
    }
  }

  void _elevateV(int t, BSplineSurface &shape)
  {
    // elevate iso-curves of constant u and then
    // reconstruct the elevated surface patch using
    // the elevated curves -- easy peasy.
    auto const n = shape.uknot.size()-shape.p-1;
    auto m = shape.vknot.size()-shape.q-1;
    // load v-const iso-curves and then elevate
    size_t i = 0;
    std::vector<BSplineCurve> curves(n,BSplineCurve());
    for (BSplineCurve &curve : curves)
    {
      curve.p = shape.q;
      curve.knot = shape.vknot;
      for (size_t j = 0; j < m; j++) curve.Q.push_back(shape.Q[shape.qid(i,j)]);
      elevate(t,curve);
      i++;
    }
    // assemble the elevated surface
    shape.q += t;
    shape.vknot = curves[0].knot;
    m = shape.vknot.size()-shape.q-1;
    shape.Q.resize(n*m);

    for (i = 0; i < n; i++)
    {
      for (size_t j = 0; j < m; j++)
      {
        shape.Q[shape.qid(i,j)] = curves[i].Q[j];
      }
    }
  }

  void elevate(size_t t, NurbsCurve &curve)
  {
    BSplineCurve b; spline_ops::weightedBSpline(curve,b);
    elevate(t,b);
    curve.p = b.p;
    curve.knot = b.knot;
    spline_ops::extractWeights(b.Q,curve.Q,curve.weights);
  }

  void elevate(size_t t, size_t direction, BSplineSurface &surf)
  {
    if (direction < 0 || direction > 1) return;
    if (direction == 0) _elevateU(t,surf);
    else                _elevateV(t,surf);
  }

  void elevate(size_t t, size_t direction, NurbsSurface &surf)
  {
    BSplineSurface b; spline_ops::weightedBSpline(surf,b);
    elevate(t,direction,b);
    surf.p = b.p;
    surf.q = b.q;
    surf.uknot = b.uknot;
    surf.vknot = b.vknot;
    spline_ops::extractWeights(b.Q,surf.Q,surf.weights);
  }
}
