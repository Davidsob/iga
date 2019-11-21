#include "DegreeReduction.h"

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

  double BezDegreeReduce(std::vector<std::vector<double>> const &bpts,
                         std::vector<std::vector<double>> &rbpts)
  {
    int const p = rbpts.size();
    int const r = (p-1)/2;
    auto alf = [p](auto const i) { return double(i)/double(p); };


    if (bpts.size() == 2)
    {
      rbpts[0] = 0.5*(bpts.back()+bpts[0]);
      return 1.0;
    }

    if ( p%2 == 0 ) // even case
    {
      rbpts[0] = bpts[0];
      rbpts.back() = bpts[p];
      for (int i = 1; i <= r; i++)
      {
        rbpts[i] = (bpts[i] - alf(i)*rbpts[i-1])/(1.0 - alf(i));
      }
      for (int i = p-2; i >= r+1; i--)
      {
        rbpts[i] = (bpts[i+1] - (1.0-alf(i+1))*rbpts[i+1])/alf(i+1);
      }
      // compute error
      return norm(bpts[r+1] - 0.5*(rbpts[r]+rbpts[r+1]));
    }
    else // odd case
    {
      rbpts[0] = bpts[0];
      rbpts.back() = bpts[p];
      for (int i = 1; i <= r-1; i++)
      {
        rbpts[i] = (bpts[i] - alf(i)*rbpts[i-1])/(1.0 - alf(i));
      }
      for (int i = p-2; i >= r+1; i--)
      {
        rbpts[i] = (bpts[i+1] - (1.0-alf(i+1))*rbpts[i+1])/alf(i+1);
      }
      // special case r
      auto const PL = (bpts[r]-alf(r)*rbpts[r-1])/(1.0-alf(r));
      auto const PR = (bpts[r+1]-(1.0-alf(r+1))*rbpts[r+1])/alf(r+1);
      rbpts[r] = 0.5*(PL + PR);

      // compute error
      return norm(PL - PR);
    }
  }

  bool reduce(BSplineCurve &curve, double tolerance)
  {
    auto const p = curve.p;
    auto const Qw = curve.Q;
    auto const U = curve.knot;
    auto const n = U.size()-p-1;

    int  const m = n+p;
    int ph = p-1; int mh = ph; int mult = p;
    int kind = ph+1; int cind = 1;
    int r = -1; int a = p; int b = p+1;

    // compute bezier degree elevation coefficeints 
    std::vector<double> const zero(curve.dim(),0);
    std::vector<double> Uh;
    std::vector<std::vector<double>> Pw;
    std::vector<std::vector<double>> bpts(p+1, zero);
    std::vector<std::vector<double>> rbpts(p, zero);
    std::vector<std::vector<double>> nextBpts(p-1, zero);
    std::vector<double> alphas(p-1,0.0);
    std::vector<double> err(m,0);

    Pw.push_back(Qw[0]);
    // compute left end of knot vector
    for (int i = 0; i <= ph; i++) Uh.push_back(U[0]);
    // initialize fist bez segment
    for (int i = 0; i <= p; i++) bpts[i] = Qw[i];

    while(b < m) // big loop through knot
    {
      // first compute knot multiplicity
      auto i = b;
      while(b < m && algo::equal(U[b],U[b+1]))
      {
        b = b+1;
      }
      mult = b-i+1;
      mh += (mult-1);
      int oldr = r; r = p-mult;
      int lbz = (oldr > 0) ? (oldr+2)/2 : 1;

      if (r > 0) // begin insert knot ub r times
      {
        auto const numer = U[b]-U[a];
        for (int k = p; k > mult; k--)
        {
          alphas[k-mult-1] = numer/(U[a+k]-U[a]);
        }

        for (int j = 1; j <= r; j++)
        {
          auto save = r-j;
          auto s = mult+j;

          for (int k = p; k >= s; k--)
          {
            double alf = alphas[k-s];
            bpts[k] = alf*bpts[k] + (1.0-alf)*bpts[k-1];
          }
          nextBpts[save] = bpts[p];
        }
      } // end insert knot

      // Degree reduce bezier segement
      auto maxErr = BezDegreeReduce(bpts,rbpts);

      err[a] += maxErr;

      if (err[a] > tolerance)
      {
        std::cout << "Curve was not degree reducible!" << std::endl;
        std::cout << "* BezDegreeReduce err = " << err[a] << std::endl;
        return false;
      }

      if (oldr > 0) // begin removing knot u = U[a]
      {
        auto first = kind; auto last = kind;

        for (int k = 0; k < oldr; k++) 
        {
          i = first; int j = last; int kj = j-kind;
          while (j-1 > k)
          {
            auto alfa = (U[a]-Uh[i-1])/(U[b]-Uh[i-1]);
            auto beta = (U[a]-Uh[j-k-1])/(U[b]-Uh[j-k-1]);
            Pw[i-1] = (Pw[i-1] - (1.0-alfa)*Pw[i-2])/alfa;
            rbpts[kj] = (rbpts[kj] - beta*rbpts[kj+1])/(1.0-beta);
            i++; j--; kj--;
          }
          // compute knot removal error bounds
          auto Br = 0.0;
          if (j-i < k)
          {
            Br = norm(Pw[i-2]-rbpts[kj+1]);
          } else {
            auto delta = (U[a]-Uh[i-1])/(U[b]-Uh[i-1]);
            auto A = delta*rbpts[kj+1]+(1.0-delta)*Pw[i-2];
            Br = norm(Pw[i-1]-A);
          }
          // update the error vector
          int K = a+oldr-k;
          int q = (2*p-k+1)/2;
          int L = K-q;

          for (int ii = L; ii <= a; ii++)
          {
            err[ii] += Br;
            if (err[ii] > tolerance)
            {
              std::cout << "Curve was not degree reducible!" << std::endl;
              std::cout << "* err = " << err[ii] << std::endl;
              return false;
            }
          }

          first--; last++;
        } // end  for (int k = 0; k < oldr; k++)

        cind--;
      } // end if (oldr > 0)

      // load knot vector and control points
      if (a != p) // load knot ua
      {
        for (i = 0; i < ph-oldr; i++)
        {
          Uh.push_back(U[a]);
          kind++;
        }
      }

      for (i = lbz; i <= ph; i++)
      {
        Pw.push_back(rbpts[i]);
        cind++;
      }
      if (b < m) // set up for next pass through loop
      {
        for (i = 0; i < r; i++)  bpts[i] = nextBpts[i];
        for (i = r; i <= p; i++) bpts[i] = Qw[b-p+i];
        a = b; b++;
      } 
      else // end knot
      {
        for (i = 0; i <= ph; i++)
        {
          Uh.push_back(U[b]);
        }
      } 
    } // end of while(b<m)

    // set new curve
    auto nr = Uh.size() - ph - 1;
    curve.p -= 1;
    curve.knot = Uh;
    curve.Q.assign(Pw.begin(), Pw.begin() + nr);
    return true;
  }

  void _reduceU(BSplineSurface &shape, double tolerance)
  {
    bool success = true;
    // reduce iso-curves of constant v and then
    // reconstruct the reduced surface patch using
    // the reduced curves -- easy peasy.
    auto const n = shape.uknot.size()-shape.p-1;
    auto const m = shape.vknot.size()-shape.q-1;
    // load v-const iso-curves and then reduce
    size_t j = 0;
    std::vector<BSplineCurve> curves(m,BSplineCurve());
    for (BSplineCurve &curve : curves)
    {
      curve.p = shape.p;
      curve.knot = shape.uknot;
      for (size_t i = 0; i < n; i++) curve.Q.push_back(shape.Q[shape.qid(i,j)]);
      success &= reduce(curve,tolerance);
      if (!success) return;
      j++;
    }

    // assemble the reduced surface
    shape.p -= 1;
    shape.uknot = curves[0].knot;
    shape.Q.clear();

    for (BSplineCurve const &curve : curves)
    {
      shape.Q.insert(shape.Q.end(),curve.Q.begin(),curve.Q.end());
    }
  }

  void _reduceV(BSplineSurface &shape, double tolerance)
  { 
    bool success = true;
    // reduce iso-curves of constant u and then
    // reconstruct the reduced surface patch using
    // the reduced curves -- easy peasy.
    auto const n = shape.uknot.size()-shape.p-1;
    auto m = shape.vknot.size()-shape.q-1;
    // load u-const iso-curves and then reduce
    size_t i = 0;
    std::vector<BSplineCurve> curves(n,BSplineCurve());
    for (BSplineCurve &curve : curves)
    {
      curve.p = shape.q;
      curve.knot = shape.vknot;
      for (size_t j = 0; j < m; j++) curve.Q.push_back(shape.Q[shape.qid(i,j)]);
      success &= reduce(curve,tolerance);
      if (!success) return;
      i++;
    }

    // assemble the reduced surface
    shape.q -= 1;
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

  template<typename Rational>
  double weightedTolerance(Rational const &shape, double const tolerance)
  {
    double wmin = 1e16;
    double Pmax = -1e16;
    for (auto const &w : shape.weights)
      wmin = std::min(wmin,std::abs(w));

    for (auto const &p : shape.Q)
      Pmax = std::max(Pmax,norm(p));

    return tolerance*wmin/(1.0+Pmax);
  }

  void reduce(NurbsCurve &curve, double tolerance)
  {
    BSplineCurve b; spline_ops::weightedBSpline(curve,b);
    auto wtol = weightedTolerance(curve,tolerance);
    reduce(b,wtol);
    curve.p = b.p;
    curve.knot = b.knot;
    spline_ops::extractWeights(b.Q,curve.Q,curve.weights);
  }

  void reduce(size_t direction, BSplineSurface &surf, double tolerance)
  {
    if (direction < 0 || direction > 1) return;
    if (direction == 0) _reduceU(surf, tolerance);
    else                _reduceV(surf, tolerance);
  }

  void reduce(size_t direction, NurbsSurface &surf, double tolerance)
  {
    BSplineSurface b; spline_ops::weightedBSpline(surf,b);
    auto wtol = weightedTolerance(surf,tolerance);
    reduce(direction,b,wtol);
    surf.p = b.p;
    surf.q = b.q;
    surf.uknot = b.uknot;
    surf.vknot = b.vknot;
    spline_ops::extractWeights(b.Q,surf.Q,surf.weights);
  }
}
