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
  inline RowVectorXd ShapeFunctions(double u, BSplineCurve const &curve)
  {
    auto n = curve.knot.size() - curve.p - 1;
    RowVectorXd shape(n); shape.setZero();

    auto const span = algo::FindSpan(u, curve.p, curve.knot);

    auto const Nu = algo::BasisFunctions(u,span,curve.p,curve.knot);

    for (int i = 0; i <= curve.p; i++)
    {
      auto a = i+span-curve.p;
      shape[a] += Nu[i];
    }

    return shape;
  }

  inline RowVectorXd ShapeFunctions(double u, NurbsCurve const &curve)
  {
    auto gen = [&curve]()
    {
      BSplineCurve b;
      b.p = curve.p;
      b.knot = curve.knot;
      b.Q = curve.Q;

      return b;
    };

    auto const shape{ShapeFunctions(u,gen())};
    auto const weights{convert::to<decltype(shape)>(curve.weights)};
    auto nrbN((shape.array()*weights.array())/shape.dot(weights));

    return  nrbN;
  }

  inline RowVectorXd ShapeFunctions(double u, double v, BSplineSurface const &surf)
  {
    auto n = surf.uknot.size() - surf.p - 1;
    auto m = surf.vknot.size() - surf.q - 1;
    RowVectorXd shape(n*m); shape.setZero();

    auto const uSpan = algo::FindSpan(u, surf.p, surf.uknot);
    auto const vSpan = algo::FindSpan(v, surf.q, surf.vknot);

    auto const Nu = algo::BasisFunctions(u,uSpan,surf.p, surf.uknot);
    auto const Nv = algo::BasisFunctions(v,vSpan,surf.q, surf.vknot);

    for (int j = 0; j <= surf.q; j++)
    {
      auto b = j+vSpan-surf.q;
      for (int i = 0; i <= surf.p; i++)
      {
        auto a = i+uSpan-surf.p;
        auto idx = surf.qid(a,b);
        shape[idx] += Nu[i]*Nv[j];
      }
    }

    return shape;
  }

  inline RowVectorXd ShapeFunctions(double u, double v, NurbsSurface const &surf)
  {
    auto gen = [&surf]()
    {
      BSplineSurface b;
      b.p = surf.p;
      b.q = surf.q;

      b.uknot = surf.uknot;
      b.vknot = surf.vknot;

      b.Q = surf.Q;

      return b;
    };

    auto const shape{ShapeFunctions(u,v,gen())};
    auto const weights{convert::to<decltype(shape)>(surf.weights)};
    auto nrbN((shape.array()*weights.array())/shape.dot(weights));

    return  nrbN;
  }

  inline RowVectorXd ShapeFunctions(double u, double v, double w, BSplineSolid const &solid)
  {
    auto n = solid.uknot.size() - solid.p - 1;
    auto m = solid.vknot.size() - solid.q - 1;
    auto o = solid.wknot.size() - solid.r - 1;
    RowVectorXd shape(n*m*o); shape.setZero();

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
  inline ShapeFunctionDerivative(double u, int order, BSplineCurve const &curve)
  {
    auto const uSpan = algo::FindSpan(u, curve.p, curve.knot);
    auto const dNu = algo::BasisFunctionDerivatives(u,uSpan,curve.p,curve.knot);

    auto const n = curve.knot.size()-curve.p-1;
    Eigen::RowVectorXd dNdU(n); dNdU.setZero();

    for (int i = 0; i <= curve.p; i++)
    {
      auto a = i+uSpan-curve.p;
      dNdU[a] += dNu[order][i];
    }

    return dNdU;
  }

  Eigen::RowVectorXd
  inline _derivativeWrtU(double u, double v, int order, BSplineSurface const &surf)
  {
    auto const uSpan = algo::FindSpan(u, surf.p, surf.uknot);
    auto const vSpan = algo::FindSpan(v, surf.q, surf.vknot);

    auto const dNu = algo::BasisFunctionDerivatives(u,uSpan,surf.p,surf.uknot);
    auto const Nv  = algo::BasisFunctions(v,vSpan,surf.q,surf.vknot);

    auto const n = surf.uknot.size()-surf.p-1;
    auto const m = surf.vknot.size()-surf.q-1;
    Eigen::RowVectorXd dNdU(n*m); dNdU.setZero();

    for (int j = 0; j <= surf.q; j++)
    {
      auto b = j+vSpan-surf.q;
      for (int i = 0; i <= surf.p; i++)
      {
        auto a = i+uSpan-surf.p;
        auto idx = surf.qid(a,b);
        dNdU[idx] += dNu[order][i]*Nv[j];
      }
    }

    return dNdU;
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

    auto const n = solid.uknot.size()-solid.p-1;
    auto const m = solid.vknot.size()-solid.q-1;
    auto const o = solid.wknot.size()-solid.r-1;
    Eigen::RowVectorXd dNdU(n*m*o); dNdU.setZero();

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
  inline _derivativeWrtV(double u, double v, int order, BSplineSurface const &surf)
  {
    auto const uSpan = algo::FindSpan(u, surf.p, surf.uknot);
    auto const vSpan = algo::FindSpan(v, surf.q, surf.vknot);

    auto const Nu  = algo::BasisFunctions(u,uSpan,surf.p,surf.uknot);
    auto const dNv = algo::BasisFunctionDerivatives(v,vSpan,surf.q,surf.vknot);

    auto const n = surf.uknot.size()-surf.p-1;
    auto const m = surf.vknot.size()-surf.q-1;
    Eigen::RowVectorXd dNdV(n*m); dNdV.setZero();

    for (int j = 0; j <= surf.q; j++)
    {
      auto b = j+vSpan-surf.q;
      for (int i = 0; i <= surf.p; i++)
      {
        auto a = i+uSpan-surf.p;
        auto idx = surf.qid(a,b);
        dNdV[idx] += Nu[i]*dNv[order][j];
      }
    }

    return dNdV;
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

    auto const n = solid.uknot.size()-solid.p-1;
    auto const m = solid.vknot.size()-solid.q-1;
    auto const o = solid.wknot.size()-solid.r-1;
    Eigen::RowVectorXd dNdV(n*m*o); dNdV.setZero();

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

    auto const n = solid.uknot.size()-solid.p-1;
    auto const m = solid.vknot.size()-solid.q-1;
    auto const o = solid.wknot.size()-solid.r-1;
    Eigen::RowVectorXd dNdW(n*m*o); dNdW.setZero();

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
  inline ShapeFunctionDerivative(double u, double v, int order, int direction, BSplineSurface const &surf)
  {
    switch (direction)
    {
      case 1:  return _derivativeWrtV(u,v,order,surf); break;
      default: return _derivativeWrtU(u,v,order,surf); break;
    }
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

  inline Eigen::RowVectorXd ShapeFunctionDerivative(double u, int order, NurbsCurve const &curve)
  {
    auto gen = [&curve]()
    {
      BSplineCurve b;
      b.p = curve.p;
      b.knot = curve.knot;
      b.Q = curve.Q;

      return b;
    };

    auto const b = gen();
    auto const shape{ShapeFunctions(u,b)};
    auto const dshape{ShapeFunctionDerivative(u,order,b)};

    auto const weights{convert::to<RowVectorXd>(curve.weights)};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const dwgt = dshape.dot(weights);

    return (weights.array()*(wgt*dshape.array() - dwgt*shape.array())/wgt2);
  }

  inline Eigen::RowVectorXd ShapeFunctionDerivative(double u, double v, int order, int direction, NurbsSurface const &surf)
  {
    auto gen = [&surf]()
    {
      BSplineSurface b;
      b.p = surf.p;
      b.q = surf.q;

      b.uknot = surf.uknot;
      b.vknot = surf.vknot;

      b.Q = surf.Q;

      return b;
    };

    auto const b = gen();
    auto const shape{ShapeFunctions(u,v,b)};
    auto const dshape{ShapeFunctionDerivative(u,v,order,direction,b)};

    auto const weights{convert::to<RowVectorXd>(surf.weights)};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const dwgt = dshape.dot(weights);

    return (weights.array()*(wgt*dshape.array() - dwgt*shape.array())/wgt2);
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

    return (weights.array()*(wgt*dshape.array() - dwgt*shape.array())/wgt2);
  }

  inline MatrixXd ShapeFunctionDerivatives(double u, NurbsCurve const &curve)
  {
    MatrixXd dN(1,curve.Q.size());
    dN.row(0) = ShapeFunctionDerivative(u,1,curve);
    return dN;
  }

  inline MatrixXd ShapeFunctionDerivatives(double u, double v, NurbsSurface const &surf)
  {
    MatrixXd dN(2,surf.Q.size());
    dN.row(0) = ShapeFunctionDerivative(u,v,1,0,surf);
    dN.row(1) = ShapeFunctionDerivative(u,v,1,1,surf);
    return dN;
  }

  inline MatrixXd ShapeFunctionDerivatives(double u, double v, double w, NurbsSolid const &solid)
  {
    MatrixXd dN(3,solid.Q.size());
    dN.row(0) = ShapeFunctionDerivative(u,v,w,1,0,solid);
    dN.row(1) = ShapeFunctionDerivative(u,v,w,1,1,solid);
    dN.row(2) = ShapeFunctionDerivative(u,v,w,1,2,solid);
    return dN;
  }

  inline Eigen::RowVectorXd parametricShapeFunction(double u)
  {
    Eigen::RowVectorXd N(2);
    N << (1.0 - u),
         (1.0 + u);

    return 0.5*N;
  }

  inline Eigen::RowVectorXd parametricShapeFunction(double u, double v)
  {
    Eigen::RowVectorXd N(4);
    N << (1.0 - u)*(1.0 - v),
         (1.0 + u)*(1.0 - v),
         (1.0 + u)*(1.0 + v),
         (1.0 - u)*(1.0 + v);

    return 0.25*N;
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

  inline Eigen::MatrixXd parametricShapeFunctionDerivatives(double u)
  {
    Eigen::RowVectorXd dN1(2);
    dN1 << -1.0, 1.0;

    Eigen::MatrixXd dN(1,2);
    dN.row(0) = dN1;
    
    return 0.5*dN;
  }

  inline Eigen::MatrixXd parametricShapeFunctionDerivatives(double u, double v)
  {
    Eigen::RowVectorXd dN1(4), dN2(4);
    dN1 << -(1.0 - v),
            (1.0 - v),
            (1.0 + v),
           -(1.0 + v);

    dN2 << -(1.0 - u),
           -(1.0 + u),
            (1.0 + u),
            (1.0 - u);

    Eigen::MatrixXd dN(2,4);
    dN.row(0) = dN1;
    dN.row(1) = dN2;
    
    return 0.25*dN;
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
    dN.row(0) = 0.125*dN1;
    dN.row(1) = 0.125*dN2;
    dN.row(2) = 0.125*dN3;
    
    return dN;
  }
}