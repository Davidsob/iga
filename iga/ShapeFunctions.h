#pragma once

#include "splines/Algorithms.h"
#include "splines/BSpline.h"
#include "splines/Nurbs.h"

#include "splines/utils/Converters.h"
#include "splines/utils/VectorOperations.h"

#include <Eigen/Dense>

using namespace Eigen;

namespace iga
{
  inline RowVectorXd ShapeFunctions(double u, double v, double w, BSplineSolid const &solid)
  {
    RowVectorXd shape(solid.Q.size()); shape.setZero();

    auto const uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    auto const vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    auto const wSpan = algo::FindSpan(w, solid.r, solid.wknot);

    auto const Nu = algo::BasisFunctions(u,uSpan,solid.p, solid.uknot);
    auto const Nv = algo::BasisFunctions(v,vSpan,solid.q, solid.vknot);
    auto const Nw = algo::BasisFunctions(w,wSpan,solid.r, solid.wknot);

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
          shape[idx] += Nu[i]*Nv[j]*Nw[k];
        }
      }
    }

    return shape;
  }

  inline RowVectorXd ShapeFunctions(double u, double v, double w, NurbsSolid const &solid)
  {
    auto gen = [&solid]()
    {
      BSplineSolid b;
      b.p = solid.p;
      b.q = solid.q;
      b.r = solid.r;

      b.uknot = solid.uknot;
      b.vknot = solid.vknot;
      b.wknot = solid.wknot;

      b.Q = solid.Q;

      return b;
    };

    auto const shape{ShapeFunctions(u,v,w,gen())};
    auto const weights{convert::to<decltype(shape)>(solid.weights)};
    auto nrbN((shape.array()*weights.array())/shape.dot(weights));

    return  nrbN;
  }

  Eigen::RowVectorXd
  inline _derivativeWrtU(double u, double v, double w, int order, BSplineSolid const &solid)
  {
    auto const uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    auto const vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    auto const wSpan = algo::FindSpan(w, solid.r, solid.wknot);

    auto const dNu = algo::BasisFunctionDerivatives(u,uSpan,solid.p,solid.uknot);
    auto const Nv  = algo::BasisFunctions(v,vSpan,solid.q,solid.vknot);
    auto const Nw  = algo::BasisFunctions(w,wSpan,solid.r,solid.wknot);

    Eigen::RowVectorXd dNdU(solid.Q.size()); dNdU.setZero();

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
          dNdU[idx] += dNu[order][i]*Nv[j]*Nw[k];
        }
      }
    }

    return dNdU;
  }

  Eigen::RowVectorXd
  inline _derivativeWrtV(double u, double v, double w, int order, BSplineSolid const &solid)
  {
    auto const uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    auto const vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    auto const wSpan = algo::FindSpan(w, solid.r, solid.wknot);

    auto const Nu  = algo::BasisFunctions(u,uSpan,solid.p,solid.uknot);
    auto const dNv = algo::BasisFunctionDerivatives(v,vSpan,solid.q,solid.vknot);
    auto const Nw  = algo::BasisFunctions(w,wSpan,solid.r,solid.wknot);

    Eigen::RowVectorXd dNdV(solid.Q.size()); dNdV.setZero();

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
          dNdV[idx] += Nu[i]*dNv[order][j]*Nw[k];
        }
      }
    }

    return dNdV;
  }

  Eigen::RowVectorXd
  inline _derivativeWrtW(double u, double v, double w, int order, BSplineSolid const &solid)
  {
    auto const uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    auto const vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    auto const wSpan = algo::FindSpan(w, solid.r, solid.wknot);

    auto const Nu  = algo::BasisFunctions(u,uSpan,solid.p,solid.uknot);
    auto const Nv  = algo::BasisFunctions(v,vSpan,solid.q,solid.vknot);
    auto const dNw = algo::BasisFunctionDerivatives(w,wSpan,solid.r,solid.wknot);

    Eigen::RowVectorXd dNdW(solid.Q.size()); dNdW.setZero();

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
          dNdW[idx] += Nu[i]*Nv[j]*dNw[order][k];
        }
      }
    }

    return dNdW;
  }

  Eigen::RowVectorXd
  inline ShapeFunctionDerivative(double u, double v, double w, int order, int direction, BSplineSolid const &solid)
  {
    switch (direction)
    {
      case 1:  return _derivativeWrtV(u,v,w,order,solid); break;
      case 2:  return _derivativeWrtW(u,v,w,order,solid); break;
      default: return _derivativeWrtU(u,v,w,order,solid); break;
    }
  }

  inline Eigen::RowVectorXd ShapeFunctionDerivative(double u, double v, double w, int order, int direction, NurbsSolid const &solid)
  {
    auto gen = [&solid]()
    {
      BSplineSolid b;
      b.p = solid.p;
      b.q = solid.q;
      b.r = solid.r;

      b.uknot = solid.uknot;
      b.vknot = solid.vknot;
      b.wknot = solid.wknot;

      b.Q = solid.Q;

      return b;
    };

    auto const b = gen();
    auto const shape{ShapeFunctions(u,v,w,b)};
    auto const dshape{ShapeFunctionDerivative(u,v,w,order,direction,b)};

    auto const weights{convert::to<RowVectorXd>(solid.weights)};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const dwgt = dshape.dot(weights);

    return (weights.array()*(wgt*dshape.array() - wgt*dwgt*shape.array())/wgt2);
  }

  inline MatrixXd ShapeFunctionDerivatives(double u, double v, double w, NurbsSolid const &solid)
  {
    MatrixXd dN(3,solid.Q.size());
    dN.row(0) = ShapeFunctionDerivative(u,v,w,1,0,solid);
    dN.row(1) = ShapeFunctionDerivative(u,v,w,1,1,solid);
    dN.row(2) = ShapeFunctionDerivative(u,v,w,1,2,solid);
    return dN;
  }

  inline Eigen::RowVectorXd parametricShapeFunction(double u, double v, double w)
  {
    Eigen::RowVectorXd N(8);
    N << (1.0 - u)*(1.0 - v)*(1.0 - w),
         (1.0 + u)*(1.0 - v)*(1.0 - w),
         (1.0 + u)*(1.0 + v)*(1.0 - w),
         (1.0 - u)*(1.0 + v)*(1.0 - w),
         (1.0 - u)*(1.0 - v)*(1.0 + w),
         (1.0 + u)*(1.0 - v)*(1.0 + w),
         (1.0 + u)*(1.0 + v)*(1.0 + w),
         (1.0 - u)*(1.0 + v)*(1.0 + w);

    return 0.125*N;
  }

  inline Eigen::MatrixXd parametricShapeFunctionDerivatives(double u, double v, double w)
  {
    Eigen::RowVectorXd dN1(8), dN2(8), dN3(8);
    dN1 << -(1.0 - v)*(1.0 - w),
            (1.0 - v)*(1.0 - w),
            (1.0 + v)*(1.0 - w),
           -(1.0 + v)*(1.0 - w),
           -(1.0 - v)*(1.0 + w),
            (1.0 - v)*(1.0 + w),
            (1.0 + v)*(1.0 + w),
           -(1.0 + v)*(1.0 + w);

    dN2 << -(1.0 - u)*(1.0 - w),
           -(1.0 + u)*(1.0 - w),
            (1.0 + u)*(1.0 - w),
            (1.0 - u)*(1.0 - w),
           -(1.0 - u)*(1.0 + w),
           -(1.0 + u)*(1.0 + w),
            (1.0 + u)*(1.0 + w),
            (1.0 - u)*(1.0 + w);

    dN3 << -(1.0 - u)*(1.0 - v),
           -(1.0 + u)*(1.0 - v),
           -(1.0 + u)*(1.0 + v),
           -(1.0 - u)*(1.0 + v),
            (1.0 - u)*(1.0 - v),
            (1.0 + u)*(1.0 - v),
            (1.0 + u)*(1.0 + v),
            (1.0 - u)*(1.0 + v);

    Eigen::MatrixXd dN(3,8);
    dN.row(0) = dN1;
    dN.row(1) = dN2;
    dN.row(2) = dN3;
    
    return 0.125*dN;
  }
}