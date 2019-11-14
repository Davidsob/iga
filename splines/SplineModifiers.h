#pragma once

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
  inline void insertKnot(double u, int r, BSplineCurve &shape)
  {
    if (r == 0) return;
    if (shape.knot.empty() || shape.Q.empty()) return;

    int const nU = shape.knot.size();
    int const nQ = shape.Q.size();
    int const  k = algo::FindSpan(u,shape.p,shape.knot);
    int const  s = algo::SpanMultiplicity(k,shape.knot);

    if ((s+r) > shape.p) r = shape.p-s; 

    // create new knot
    decltype(shape.knot) U(nU+r,0.0);
    for (int i = 0; i <= k; i++) U[i] = shape.knot[i];
    for (int i = 1; i <= r; i++) U[i+k] = u;
    for (int i = k+1; i < nU; i++) U[i+r] = shape.knot[i];

    // save unaltered control points
    std::vector<double> const zero(shape.dim(),0);
    std::vector<std::vector<double>> tmp(shape.p+1,zero);
    decltype(shape.Q) P(nQ+r,zero);

    for (int i = 0; i <= k-shape.p; i++) P[i] = shape.Q[i];
    for (int i = k; i < nQ; i++) P[i+r] = shape.Q[i];
    for (int i = 0; i <= shape.p-s; i++) tmp[i] = shape.Q[k-shape.p+i];

    int L(-1);
    for (int j = 1; j <= r; j++)
    {
      L = k-shape.p+j;
      for (int i = 0; i <= shape.p-j-s; i++)
      {
        auto alf = (u-shape.knot[L+i])/(shape.knot[i+k+1] - shape.knot[L+i]);
        tmp[i] = alf*tmp[i+1] + (1.0-alf)*tmp[i];
      }

      P[L] = tmp[0];
      P[k+r-j] = tmp[shape.p-j];
    }

    // load the remaining control points
    for (int i = L+1; i < k; i++) P[i] = tmp[i-L];

    // update the shape
    shape.knot = U;
    shape.Q = P;
  }

  inline void insertKnot(double uv, int r, int direction, BSplineSurface &shape)
  {
    if (r == 0 || shape.Q.empty()) return;

    int const nQ = shape.Q.size();
    int const  N = shape.p+1;
    int const  M = shape.q+1;

    if (direction == 0)
    {
      if (shape.uknot.empty()) return;

      int const nU = shape.uknot.size();
      int const  k = algo::FindSpan(uv,shape.p,shape.uknot);
      int const  s = algo::SpanMultiplicity(k,shape.uknot);

      if ((s+r) > shape.p) r = shape.p-s; 

      // create new uknot
      decltype(shape.uknot) U(nU+r,0.0);

      for (int i = 0; i <= k; i++) U[i] = shape.uknot[i];
      for (int i = 1; i <= r; i++) U[i+k] = uv;
      for (int i = k+1; i < nU; i++) U[i+r] = shape.uknot[i];

      // alpha table
      std::vector<std::vector<double>> alf(r+1,std::vector<double>(r+1,0));
      int L(-1);
      for (int j = 1; j <= r; j++)
      {
        L = k-shape.p+j;
        for (int i = 0; i <= shape.p-j-s; i++)
        {
          alf[i][j] = (uv-U[L+i])/(U[i+k+1] - U[L+i]);
        }
      }

      std::vector<double> zero(shape.dim(),0);
      decltype(shape.Q) P(nQ+M*r,zero);
      int const Un = U.size()-shape.p-1;
      auto const idx = [Un](auto i, auto j) { return i + j*Un; };
      for (int row = 0; row < M; row++)
      {
        // save unaltered control points
        decltype(shape.Q) tmp(shape.p+1,zero);
        for (int i = 0; i <= k-shape.p; i++) P[idx(i,row)] = shape.Q[shape.qid(i,row)];
        for (int i = k; i < N; i++) P[idx(i+r,row)] = shape.Q[shape.qid(i,row)];
        for (int i = 0; i <= shape.p-s; i++) tmp[i] = shape.Q[shape.qid(k-shape.p+i,row)];

        for (int j = 1; j <= r; j++)
        {
          L = k-shape.p+j;
          for (int i = 0; i <= shape.p-j-s; i++)
          {
            auto fact = alf[i][j];
            tmp[i] = fact*tmp[i+1] + (1.0-fact)*tmp[i];
          }

          P[idx(L,row)] = tmp[0];
          P[idx(k+r-j,row)] = tmp[shape.p-j];
        }
        // load the remaining control points
        for (int i = L+1; i < k; i++) P[idx(i,row)] = tmp[i-L];
      }

      // update the shape
      shape.uknot = U;
      shape.Q = P;
    } else if (direction == 1)
    {

      if (shape.vknot.empty()) return;

      int const nV = shape.vknot.size();
      int const  k = algo::FindSpan(uv,shape.q,shape.vknot);
      int const  s = algo::SpanMultiplicity(k,shape.vknot);

      if ((s+r) > shape.q) r = shape.q-s; 

      // create new uknot
      decltype(shape.vknot) V(nV+r,0.0);

      for (int i = 0; i <= k; i++) V[i] = shape.vknot[i];
      for (int i = 1; i <= r; i++) V[i+k] = uv;
      for (int i = k+1; i < nV; i++) V[i+r] = shape.vknot[i];

      // alpha table
      int L(-1);
      std::vector<std::vector<double>> alf(r+1,std::vector<double>(r+1,0));
      std::cout << V << std::endl;
      for (int j = 1; j <= r; j++)
      {
        L = k-shape.q+j;
        for (int i = 0; i <= r; i++)
        {
          alf[i][j] = (uv-V[L+i])/(V[i+k+1] - V[L+i]);
        }
      }
      std::cout << alf << std::endl;
      
      std::vector<double> zero(shape.dim(),0);
      decltype(shape.Q) P(nQ+N*r,zero);
      int const Un = shape.uknot.size()-shape.p-1;
      auto const idx = [Un](auto i, auto j) { return i + j*Un; };
      for (int col = 0; col < N; col++)
      {
        // save unaltered control points
        decltype(shape.Q) tmp(shape.q+1,zero);
        for (int i = 0; i <= k-shape.q; i++) P[idx(col,i)] = shape.Q[shape.qid(col,i)];
        for (int i = k; i < M; i++) P[idx(col,i+r)] = shape.Q[shape.qid(col,i)];
        for (int i = 0; i <= shape.q-s; i++) tmp[i] = shape.Q[shape.qid(col,k-shape.q+i)];

        for (int j = 1; j <= r; j++)
        {
          L = k-shape.q+j;
          for (int i = 0; i <= shape.q-j-s; i++)
          {
            auto fact = alf[i][j];
            tmp[i] = fact*tmp[i+1] + (1.0-fact)*tmp[i];
          }

          P[idx(col,L)] = tmp[0];
          P[idx(col,k+r-j)] = tmp[shape.q-j];
        }
        std::cout << alf << std::endl;
        std::cout << P << std::endl;
        // // load the remaining control points
        // for (int i = L+1; i < k; i++) P[idx(col,i)] = tmp[i-L];
      }

      // update the shape
      shape.vknot = V;
      shape.Q = P;
    }
  }

  inline void insertKnot(double u, size_t rep, NurbsCurve &curve)
  {
    BSplineCurve b; spline_ops::weightedBSpline(curve,b);
    insertKnot(u,rep,b);
    curve.knot = b.knot;
    spline_ops::extractWeights(b.Q,curve.Q,curve.weights);
  }

  inline void insertKnot(double uv, size_t rep, size_t direction, NurbsSurface &surface)
  {
    BSplineSurface b; spline_ops::weightedBSpline(surface,b);
    insertKnot(uv,rep,direction,b);
    surface.uknot = b.uknot;
    surface.vknot = b.vknot;
    spline_ops::extractWeights(b.Q,surface.Q,surface.weights);
  }
}