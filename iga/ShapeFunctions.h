#pragma once

#include "splines/Algorithms.h"
#include "splines/BSpline.h"
#include "splines/Nurbs.h"

#include "splines/utils/Converters.h"
#include "splines/utils/VectorOperations.h"

#include <Eigen/Dense>

namespace iga
{
  template<typename Curve>
  inline std::vector<int> ActiveControlPoints(double u, Curve const &curve)
  {
    auto const span = algo::FindSpan(u, curve.p, curve.knot);

    std::vector<int> active;
    for (int i = 0; i <= curve.p; i++)
    {
      int a = i+span-curve.p;
      active.push_back(a);
    }
    return active;
  }

  template<typename Surface>
  inline std::vector<int> ActiveControlPoints(double u, double v, Surface const &surf)
  {
    auto const uSpan = algo::FindSpan(u, surf.p, surf.uknot);
    auto const vSpan = algo::FindSpan(v, surf.q, surf.vknot);

    std::vector<int> active;
    for (int j = 0; j <= surf.q; j++)
    {
      auto b = j+vSpan-surf.q;
      for (int i = 0; i <= surf.p; i++)
      {
        auto a = i+uSpan-surf.p;
        auto idx = surf.qid(a,b);
        active.push_back(idx);
      }
    }
    return active;
  }

  template<typename Solid>
  inline std::vector<int> ActiveControlPoints(double u, double v, double w, Solid const &solid)
  {
    auto const uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    auto const vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    auto const wSpan = algo::FindSpan(w, solid.r, solid.wknot);

    std::vector<int> active;
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
          active.push_back(idx);
        }
      }
    }

    return active;
  }

  inline Eigen::RowVectorXd CompactShapeFunctions(double u, BSplineCurve const &curve)
  {
    auto const span = algo::FindSpan(u, curve.p, curve.knot);
    auto const Nu = algo::BasisFunctions(u,span,curve.p,curve.knot);
    return convert::to<Eigen::RowVectorXd>(Nu);
  }

  inline Eigen::RowVectorXd ShapeFunctions(double u, BSplineCurve const &curve)
  {
    auto n = curve.knot.size() - curve.p - 1;
    Eigen::RowVectorXd shape(n); shape.setZero();

    auto const span = algo::FindSpan(u, curve.p, curve.knot);

    auto const Nu = algo::BasisFunctions(u,span,curve.p,curve.knot);

    for (int i = 0; i <= curve.p; i++)
    {
      auto a = i+span-curve.p;
      shape[a] += Nu[i];
    }

    return shape;
  }

  inline Eigen::RowVectorXd ShapeFunctions(double u, NurbsCurve const &curve)
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

  inline Eigen::RowVectorXd CompactShapeFunctions(double u, NurbsCurve const &curve)
  {
    auto gen = [&curve]()
    {
      BSplineCurve b;
      b.p = curve.p;
      b.knot = curve.knot;
      b.Q = curve.Q;

      return b;
    };

    auto const bcurve = gen();
    auto const active{ActiveControlPoints(u,bcurve)};
    auto const shape{CompactShapeFunctions(u,bcurve)};
    auto const weights{convert::to<decltype(shape)>(subVector(curve.weights,active))};
    auto nrbN((shape.array()*weights.array())/shape.dot(weights));

    return  nrbN;
  }

  inline Eigen::RowVectorXd CompactShapeFunctions(double u, double v, BSplineSurface const &surf)
  {
    std::vector<double> shape;

    auto const uSpan = algo::FindSpan(u, surf.p, surf.uknot);
    auto const vSpan = algo::FindSpan(v, surf.q, surf.vknot);

    auto const Nu = algo::BasisFunctions(u,uSpan,surf.p, surf.uknot);
    auto const Nv = algo::BasisFunctions(v,vSpan,surf.q, surf.vknot);

    for (int j = 0; j <= surf.q; j++)
    {
      for (int i = 0; i <= surf.p; i++)
      {
        shape.push_back(Nu[i]*Nv[j]);
      }
    }

    return convert::to<Eigen::RowVectorXd>(shape);
  }

  inline Eigen::RowVectorXd ShapeFunctions(double u, double v, BSplineSurface const &surf)
  {
    auto n = surf.uknot.size() - surf.p - 1;
    auto m = surf.vknot.size() - surf.q - 1;
    Eigen::RowVectorXd shape(n*m); shape.setZero();

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

  inline Eigen::RowVectorXd CompactShapeFunctions(double u, double v, NurbsSurface const &surf)
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

    auto const bsurf = gen();
    auto const active{ActiveControlPoints(u,v,bsurf)};
    auto const shape{CompactShapeFunctions(u,v,bsurf)};
    auto const weights{convert::to<decltype(shape)>(subVector(surf.weights,active))};
    auto nrbN((shape.array()*weights.array())/shape.dot(weights));

    return  nrbN;
  }

  inline Eigen::RowVectorXd ShapeFunctions(double u, double v, NurbsSurface const &surf)
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

  inline Eigen::RowVectorXd CompactShapeFunctions(double u, double v, double w, BSplineSolid const &solid)
  {
    std::vector<double> shape;

    auto const uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    auto const vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    auto const wSpan = algo::FindSpan(w, solid.r, solid.wknot);

    auto const Nu = algo::BasisFunctions(u,uSpan,solid.p, solid.uknot);
    auto const Nv = algo::BasisFunctions(v,vSpan,solid.q, solid.vknot);
    auto const Nw = algo::BasisFunctions(w,wSpan,solid.r, solid.wknot);

    for (int k = 0; k <= solid.r; k++)
    {
      for (int j = 0; j <= solid.q; j++)
      {
        for (int i = 0; i <= solid.p; i++)
        {
          shape.push_back(Nu[i]*Nv[j]*Nw[k]);
        }
      }
    }

    return convert::to<Eigen::RowVectorXd>(shape);
  }

  inline Eigen::RowVectorXd ShapeFunctions(double u, double v, double w, BSplineSolid const &solid)
  {
    auto n = solid.uknot.size() - solid.p - 1;
    auto m = solid.vknot.size() - solid.q - 1;
    auto o = solid.wknot.size() - solid.r - 1;
    Eigen::RowVectorXd shape(n*m*o); shape.setZero();

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

  inline Eigen::RowVectorXd CompactShapeFunctions(double u, double v, double w, NurbsSolid const &solid)
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

    auto const bsolid = gen();
    auto const active{ActiveControlPoints(u,v,w,bsolid)};
    auto const shape{ShapeFunctions(u,v,w,bsolid)};
    auto const weights{convert::to<decltype(shape)>(subVector(solid.weights,active))};
    auto nrbN((shape.array()*weights.array())/shape.dot(weights));

    return  nrbN;
  }

  inline Eigen::RowVectorXd ShapeFunctions(double u, double v, double w, NurbsSolid const &solid)
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

  inline Eigen::RowVectorXd
  CompactShapeFunctionDerivative(double u, int order, BSplineCurve const &curve)
  {
    auto const uSpan = algo::FindSpan(u, curve.p, curve.knot);
    auto const dNu = algo::BasisFunctionDerivatives(u,uSpan,curve.p,curve.knot);

    std::vector<double> dNdU;
    for (int i = 0; i <= curve.p; i++)
    {
      dNdU.push_back(dNu[order][i]);
    }

    return convert::to<Eigen::RowVectorXd>(dNdU);
  }

  inline Eigen::RowVectorXd
  ShapeFunctionDerivative(double u, int order, BSplineCurve const &curve)
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

  inline Eigen::RowVectorXd
  _CompactDerivativeWrtU(double u, double v, int order, BSplineSurface const &surf)
  {
    auto const uSpan = algo::FindSpan(u, surf.p, surf.uknot);
    auto const vSpan = algo::FindSpan(v, surf.q, surf.vknot);

    auto const dNu = algo::BasisFunctionDerivatives(u,uSpan,surf.p,surf.uknot);
    auto const Nv  = algo::BasisFunctions(v,vSpan,surf.q,surf.vknot);

    std::vector<double> dNdU;

    for (int j = 0; j <= surf.q; j++)
    {
      for (int i = 0; i <= surf.p; i++)
      {
        dNdU.push_back(dNu[order][i]*Nv[j]);
      }
    }

    return convert::to<Eigen::RowVectorXd>(dNdU);
  }

  inline Eigen::RowVectorXd
  _derivativeWrtU(double u, double v, int order, BSplineSurface const &surf)
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

  inline Eigen::RowVectorXd
  _CompactDerivativeWrtU(double u, double v, double w, int order, BSplineSolid const &solid)
  {
    auto const uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    auto const vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    auto const wSpan = algo::FindSpan(w, solid.r, solid.wknot);

    auto const dNu = algo::BasisFunctionDerivatives(u,uSpan,solid.p,solid.uknot);
    auto const Nv  = algo::BasisFunctions(v,vSpan,solid.q,solid.vknot);
    auto const Nw  = algo::BasisFunctions(w,wSpan,solid.r,solid.wknot);

    std::vector<double> dNdU;
    for (int k = 0; k <= solid.r; k++)
    {
      for (int j = 0; j <= solid.q; j++)
      {
        for (int i = 0; i <= solid.p; i++)
        {
          dNdU.push_back(dNu[order][i]*Nv[j]*Nw[k]);
        }
      }
    }

    return convert::to<Eigen::RowVectorXd>(dNdU);
  }

  inline Eigen::RowVectorXd
  _derivativeWrtU(double u, double v, double w, int order, BSplineSolid const &solid)
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

  inline Eigen::RowVectorXd
  CompactShapeFunctionMixedDerivative(double u, double v, BSplineSurface const &surf)
  {
    auto const uSpan = algo::FindSpan(u, surf.p, surf.uknot);
    auto const vSpan = algo::FindSpan(v, surf.q, surf.vknot);

    auto const dNu = algo::BasisFunctionDerivatives(u,uSpan,surf.p,surf.uknot);
    auto const dNv = algo::BasisFunctionDerivatives(v,vSpan,surf.q,surf.vknot);

    std::vector<double> dNdUV;
    for (int j = 0; j <= surf.q; j++)
    {
      for (int i = 0; i <= surf.p; i++)
      {
        dNdUV.push_back(dNu[1][i]*dNv[1][j]);
      }
    }

    return convert::to<Eigen::RowVectorXd>(dNdUV);
  }

  inline Eigen::RowVectorXd
  ShapeFunctionMixedDerivative(double u, double v, BSplineSurface const &surf)
  {
    auto const uSpan = algo::FindSpan(u, surf.p, surf.uknot);
    auto const vSpan = algo::FindSpan(v, surf.q, surf.vknot);

    auto const dNu = algo::BasisFunctionDerivatives(u,uSpan,surf.p,surf.uknot);
    auto const dNv = algo::BasisFunctionDerivatives(v,vSpan,surf.q,surf.vknot);

    auto const n = surf.uknot.size()-surf.p-1;
    auto const m = surf.vknot.size()-surf.q-1;
    Eigen::RowVectorXd dNdUV(n*m); dNdUV.setZero();

    for (int j = 0; j <= surf.q; j++)
    {
      auto b = j+vSpan-surf.q;
      for (int i = 0; i <= surf.p; i++)
      {
        auto a = i+uSpan-surf.p;
        auto idx = surf.qid(a,b);
        dNdUV[idx] += dNu[1][i]*dNv[1][j];
      }
    }

    return dNdUV;
  }

  inline Eigen::RowVectorXd
  _CompactDerivativeWrtV(double u, double v, int order, BSplineSurface const &surf)
  {
    auto const uSpan = algo::FindSpan(u, surf.p, surf.uknot);
    auto const vSpan = algo::FindSpan(v, surf.q, surf.vknot);

    auto const Nu  = algo::BasisFunctions(u,uSpan,surf.p,surf.uknot);
    auto const dNv = algo::BasisFunctionDerivatives(v,vSpan,surf.q,surf.vknot);

    std::vector<double> dNdV;
    for (int j = 0; j <= surf.q; j++)
    {
      for (int i = 0; i <= surf.p; i++)
      {
        dNdV.push_back(Nu[i]*dNv[order][j]);
      }
    }

    return convert::to<Eigen::RowVectorXd>(dNdV);
  }

  inline Eigen::RowVectorXd
  _derivativeWrtV(double u, double v, int order, BSplineSurface const &surf)
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

  inline Eigen::RowVectorXd
  _CompactDerivativeWrtV(double u, double v, double w, int order, BSplineSolid const &solid)
  {
    auto const uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    auto const vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    auto const wSpan = algo::FindSpan(w, solid.r, solid.wknot);

    auto const Nu  = algo::BasisFunctions(u,uSpan,solid.p,solid.uknot);
    auto const dNv = algo::BasisFunctionDerivatives(v,vSpan,solid.q,solid.vknot);
    auto const Nw  = algo::BasisFunctions(w,wSpan,solid.r,solid.wknot);

    std::vector<double> dNdV;
    for (int k = 0; k <= solid.r; k++)
    {
      for (int j = 0; j <= solid.q; j++)
      {
        for (int i = 0; i <= solid.p; i++)
        {
          dNdV.push_back(Nu[i]*dNv[order][j]*Nw[k]);
        }
      }
    }

    return convert::to<Eigen::RowVectorXd>(dNdV);
  }

  inline Eigen::RowVectorXd
  _derivativeWrtV(double u, double v, double w, int order, BSplineSolid const &solid)
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

  inline Eigen::RowVectorXd
  _CompactDerivativeWrtW(double u, double v, double w, int order, BSplineSolid const &solid)
  {
    auto const uSpan = algo::FindSpan(u, solid.p, solid.uknot);
    auto const vSpan = algo::FindSpan(v, solid.q, solid.vknot);
    auto const wSpan = algo::FindSpan(w, solid.r, solid.wknot);

    auto const Nu  = algo::BasisFunctions(u,uSpan,solid.p,solid.uknot);
    auto const Nv  = algo::BasisFunctions(v,vSpan,solid.q,solid.vknot);
    auto const dNw = algo::BasisFunctionDerivatives(w,wSpan,solid.r,solid.wknot);

    std::vector<double> dNdW;
    for (int k = 0; k <= solid.r; k++)
    {
      for (int j = 0; j <= solid.q; j++)
      {
        for (int i = 0; i <= solid.p; i++)
        {
          dNdW.push_back(Nu[i]*Nv[j]*dNw[order][k]);
        }
      }
    }

    return convert::to<Eigen::RowVectorXd>(dNdW);
  }

  inline Eigen::RowVectorXd
  _derivativeWrtW(double u, double v, double w, int order, BSplineSolid const &solid)
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

  inline Eigen::RowVectorXd
  CompactShapeFunctionDerivative(double u, double v, int order, int direction, BSplineSurface const &surf)
  {
    switch (direction)
    {
      case 1:  return _CompactDerivativeWrtV(u,v,order,surf); break;
      default: return _CompactDerivativeWrtU(u,v,order,surf); break;
    }
  }

  inline Eigen::RowVectorXd
  ShapeFunctionDerivative(double u, double v, int order, int direction, BSplineSurface const &surf)
  {
    switch (direction)
    {
      case 1:  return _derivativeWrtV(u,v,order,surf); break;
      default: return _derivativeWrtU(u,v,order,surf); break;
    }
  }

  inline Eigen::RowVectorXd
  CompactShapeFunctionDerivative(double u, double v, double w, int order, int direction, BSplineSolid const &solid)
  {
    switch (direction)
    {
      case 1:  return _CompactDerivativeWrtV(u,v,w,order,solid); break;
      case 2:  return _CompactDerivativeWrtW(u,v,w,order,solid); break;
      default: return _CompactDerivativeWrtU(u,v,w,order,solid); break;
    }
  }

  inline Eigen::RowVectorXd
  ShapeFunctionDerivative(double u, double v, double w, int order, int direction, BSplineSolid const &solid)
  {
    switch (direction)
    {
      case 1:  return _derivativeWrtV(u,v,w,order,solid); break;
      case 2:  return _derivativeWrtW(u,v,w,order,solid); break;
      default: return _derivativeWrtU(u,v,w,order,solid); break;
    }
  }

  inline Eigen::RowVectorXd
  CompactShapeFunctionDerivative(double u, NurbsCurve const &curve)
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
    auto const active{ActiveControlPoints(u,b)};
    auto const shape{CompactShapeFunctions(u,b)};
    auto const dshape{CompactShapeFunctionDerivative(u,1,b)};

    auto const weights{convert::to<Eigen::RowVectorXd>(subVector(curve.weights,active))};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const dwgt = dshape.dot(weights);

    return (weights.array()*(wgt*dshape.array() - dwgt*shape.array())/wgt2);
  }

  inline Eigen::RowVectorXd
  ShapeFunctionDerivative(double u, NurbsCurve const &curve)
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
    auto const dshape{ShapeFunctionDerivative(u,1,b)};

    auto const weights{convert::to<Eigen::RowVectorXd>(curve.weights)};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const dwgt = dshape.dot(weights);

    return (weights.array()*(wgt*dshape.array() - dwgt*shape.array())/wgt2);
  }

  inline Eigen::RowVectorXd
  CompactShapeFunctionDerivative2(double u, NurbsCurve const &curve)
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
    auto const active{ActiveControlPoints(u,b)};
    auto const shape{CompactShapeFunctions(u,b)};
    auto const dshape{CompactShapeFunctionDerivative(u,1,b)};
    auto const dshape2{CompactShapeFunctionDerivative(u,2,b)};

    auto const weights{convert::to<Eigen::RowVectorXd>(subVector(curve.weights,active))};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const wgt3 = wgt2*wgt;
    auto const dwgt = dshape.dot(weights);
    auto const dwgt2 = dshape2.dot(weights);

    auto A = (weights.array()*(wgt*dshape2.array() - dwgt2*shape.array())/wgt2);
    auto B = 2.0*dwgt*(weights.array()*(wgt*dshape.array() - dwgt*shape.array()))/wgt3;
    return A - B;
  }

  inline Eigen::RowVectorXd
  ShapeFunctionDerivative2(double u, NurbsCurve const &curve)
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
    auto const dshape{ShapeFunctionDerivative(u,1,b)};
    auto const dshape2{ShapeFunctionDerivative(u,2,b)};

    auto const weights{convert::to<Eigen::RowVectorXd>(curve.weights)};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const wgt3 = wgt2*wgt;
    auto const dwgt = dshape.dot(weights);
    auto const dwgt2 = dshape2.dot(weights);

    auto A = (weights.array()*(wgt*dshape2.array() - dwgt2*shape.array())/wgt2);
    auto B = 2.0*dwgt*(weights.array()*(wgt*dshape.array() - dwgt*shape.array()))/wgt3;
    return A - B;
  }

  inline Eigen::RowVectorXd
  CompactShapeFunctionDerivative(double u, double v, int direction, NurbsSurface const &surf)
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
    auto const active{ActiveControlPoints(u,v,b)};
    auto const shape{CompactShapeFunctions(u,v,b)};
    auto const dshape{CompactShapeFunctionDerivative(u,v,1,direction,b)};

    auto const weights{convert::to<Eigen::RowVectorXd>(subVector(surf.weights,active))};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const dwgt = dshape.dot(weights);

    return (weights.array()*(wgt*dshape.array() - dwgt*shape.array())/wgt2);
  }

  inline Eigen::RowVectorXd
  ShapeFunctionDerivative(double u, double v, int direction, NurbsSurface const &surf)
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
    auto const dshape{ShapeFunctionDerivative(u,v,1,direction,b)};

    auto const weights{convert::to<Eigen::RowVectorXd>(surf.weights)};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const dwgt = dshape.dot(weights);

    return (weights.array()*(wgt*dshape.array() - dwgt*shape.array())/wgt2);
  }

  inline Eigen::RowVectorXd
  CompactShapeFunctionDerivative2(double u, double v, int direction, NurbsSurface const &surf)
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
    auto const active{ActiveControlPoints(u,v,b)};
    auto const shape  {CompactShapeFunctions(u,v,b)};
    auto const dshape {CompactShapeFunctionDerivative(u,v,1,direction,b)};
    auto const dshape2{CompactShapeFunctionDerivative(u,v,2,direction,b)};

    auto const weights{convert::to<Eigen::RowVectorXd>(subVector(surf.weights,active))};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const wgt3 = wgt2*wgt;
    auto const dwgt = dshape.dot(weights);
    auto const dwgt2 = dshape2.dot(weights);

    auto A = (weights.array()*(wgt*dshape2.array() - dwgt2*shape.array())/wgt2);
    auto B = 2.0*dwgt*(weights.array()*(wgt*dshape.array() - dwgt*shape.array()))/wgt3;
    return A - B;
  }

  inline Eigen::RowVectorXd
  ShapeFunctionDerivative2(double u, double v, int direction, NurbsSurface const &surf)
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
    auto const dshape{ShapeFunctionDerivative(u,v,1,direction,b)};
    auto const dshape2{ShapeFunctionDerivative(u,v,2,direction,b)};

    auto const weights{convert::to<Eigen::RowVectorXd>(surf.weights)};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const wgt3 = wgt2*wgt;
    auto const dwgt = dshape.dot(weights);
    auto const dwgt2 = dshape2.dot(weights);

    auto A = (weights.array()*(wgt*dshape2.array() - dwgt2*shape.array())/wgt2);
    auto B = 2.0*dwgt*(weights.array()*(wgt*dshape.array() - dwgt*shape.array()))/wgt3;
    return A - B;
  }

  inline Eigen::RowVectorXd
  CompactShapeFunctionMixedDerivative(double u, double v, NurbsSurface const &surf)
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
    auto const active{ActiveControlPoints(u,v,b)};
    auto const shape   {CompactShapeFunctions(u,v,b)};
    auto const dshapeU {CompactShapeFunctionDerivative(u,v,1,0,b)};
    auto const dshapeV {CompactShapeFunctionDerivative(u,v,1,1,b)};
    auto const dshapeUV{CompactShapeFunctionMixedDerivative(u,v,b)};

    auto const weights{convert::to<Eigen::RowVectorXd>(subVector(surf.weights,active))};
    auto const wgt    = shape.dot(weights);
    auto const wgt2   = std::pow(wgt,2);
    auto const wgt3   = wgt2*wgt;
    auto const dwgtU  = dshapeU.dot(weights);
    auto const dwgtV  = dshapeV.dot(weights);
    auto const dwgtUV = dshapeUV.dot(weights);

    Eigen::RowVectorXd Nuv =
      weights.array()*(wgt2*dshapeUV.array()
                     - wgt*(dshapeU.array()*dwgtV + dshapeV.array()*dwgtU + shape.array()*dwgtUV)
                     + 2.0*shape.array()*dwgtU*dwgtV)/wgt3;
    return Nuv;
  }

  inline Eigen::RowVectorXd
  ShapeFunctionMixedDerivative(double u, double v, NurbsSurface const &surf)
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
    auto const dshapeU{ShapeFunctionDerivative(u,v,1,0,b)};
    auto const dshapeV{ShapeFunctionDerivative(u,v,1,1,b)};
    auto const dshapeUV{ShapeFunctionMixedDerivative(u,v,b)};

    auto const weights{convert::to<Eigen::RowVectorXd>(surf.weights)};
    auto const wgt    = shape.dot(weights);
    auto const wgt2   = std::pow(wgt,2);
    auto const wgt3   = wgt2*wgt;
    auto const dwgtU  = dshapeU.dot(weights);
    auto const dwgtV  = dshapeV.dot(weights);
    auto const dwgtUV = dshapeUV.dot(weights);

    Eigen::RowVectorXd Nuv =
      weights.array()*(wgt2*dshapeUV.array()
                     - wgt*(dshapeU.array()*dwgtV + dshapeV.array()*dwgtU + shape.array()*dwgtUV)
                     + 2.0*shape.array()*dwgtU*dwgtV)/wgt3;
    return Nuv;
  }

  inline Eigen::RowVectorXd
  CompactShapeFunctionDerivative(double u, double v, double w, int direction, NurbsSolid const &solid)
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
    auto const active{ActiveControlPoints(u,v,w,b)};
    auto const shape {CompactShapeFunctions(u,v,w,b)};
    auto const dshape{CompactShapeFunctionDerivative(u,v,w,1,direction,b)};

    auto const weights{convert::to<Eigen::RowVectorXd>(subVector(solid.weights,active))};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const dwgt = dshape.dot(weights);

    return (weights.array()*(wgt*dshape.array() - dwgt*shape.array())/wgt2);
  }

  inline Eigen::RowVectorXd
  ShapeFunctionDerivative(double u, double v, double w, int direction, NurbsSolid const &solid)
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
    auto const dshape{ShapeFunctionDerivative(u,v,w,1,direction,b)};

    auto const weights{convert::to<Eigen::RowVectorXd>(solid.weights)};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const dwgt = dshape.dot(weights);

    return (weights.array()*(wgt*dshape.array() - dwgt*shape.array())/wgt2);
  }

  inline Eigen::RowVectorXd
  CompactShapeFunctionDerivative2(double u, double v, double w, int direction, NurbsSolid const &solid)
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
    auto const active{ActiveControlPoints(u,v,w,b)};
    auto const shape  {CompactShapeFunctions(u,v,w,b)};
    auto const dshape {CompactShapeFunctionDerivative(u,v,w,1,direction,b)};
    auto const dshape2{CompactShapeFunctionDerivative(u,v,w,2,direction,b)};

    auto const weights{convert::to<Eigen::RowVectorXd>(subVector(solid.weights,active))};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const wgt3 = wgt2*wgt;
    auto const dwgt = dshape.dot(weights);
    auto const dwgt2 = dshape2.dot(weights);

    auto A = (weights.array()*(wgt*dshape2.array() - dwgt2*shape.array())/wgt2);
    auto B = 2.0*dwgt*(weights.array()*(wgt*dshape.array() - dwgt*shape.array()))/wgt3;
    return A - B;
  }

  inline Eigen::RowVectorXd
  ShapeFunctionDerivative2(double u, double v, double w, int direction, NurbsSolid const &solid)
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
    auto const dshape{ShapeFunctionDerivative(u,v,w,1,direction,b)};
    auto const dshape2{ShapeFunctionDerivative(u,v,w,2,direction,b)};

    auto const weights{convert::to<Eigen::RowVectorXd>(solid.weights)};
    auto const wgt  = shape.dot(weights);
    auto const wgt2 = std::pow(wgt,2);
    auto const wgt3 = wgt2*wgt;
    auto const dwgt = dshape.dot(weights);
    auto const dwgt2 = dshape2.dot(weights);

    auto A = (weights.array()*(wgt*dshape2.array() - dwgt2*shape.array())/wgt2);
    auto B = 2.0*dwgt*(weights.array()*(wgt*dshape.array() - dwgt*shape.array()))/wgt3;
    return A - B;
  }

  inline Eigen::MatrixXd CompactShapeFunctionDerivatives(double u, NurbsCurve const &curve)
  {
    Eigen::RowVectorXd const dNd1 = CompactShapeFunctionDerivative(u,curve);
    Eigen::MatrixXd dN(1,dNd1.size()); dN << dNd1;
    return dN;
  }

  inline Eigen::MatrixXd ShapeFunctionDerivatives(double u, NurbsCurve const &curve)
  {
    Eigen::MatrixXd dN(1,curve.Q.size());
    dN.row(0) = ShapeFunctionDerivative(u,curve);
    return dN;
  }

  inline Eigen::MatrixXd CompactShapeFunctionDerivatives(double u, double v, NurbsSurface const &surf)
  {
    Eigen::RowVectorXd const dNd1 = CompactShapeFunctionDerivative(u,v,0,surf);
    Eigen::RowVectorXd const dNd2 = CompactShapeFunctionDerivative(u,v,1,surf);
    Eigen::MatrixXd dN(2,dNd1.size());
    dN.row(0) = dNd1;
    dN.row(1) = dNd2;
    return dN;
  }

  inline Eigen::MatrixXd ShapeFunctionDerivatives(double u, double v, NurbsSurface const &surf)
  {
    Eigen::MatrixXd dN(2,surf.Q.size());
    dN.row(0) = ShapeFunctionDerivative(u,v,0,surf);
    dN.row(1) = ShapeFunctionDerivative(u,v,1,surf);
    return dN;
  }

  inline Eigen::MatrixXd CompactShapeFunctionDerivatives(double u, double v, double w, NurbsSolid const &solid)
  {
    Eigen::RowVectorXd const dNd1 = CompactShapeFunctionDerivative(u,v,w,0,solid);
    Eigen::RowVectorXd const dNd2 = CompactShapeFunctionDerivative(u,v,w,1,solid);
    Eigen::RowVectorXd const dNd3 = CompactShapeFunctionDerivative(u,v,w,2,solid);
    Eigen::MatrixXd dN(3,dNd1.size());
    dN.row(0) = dNd1;
    dN.row(1) = dNd2;
    dN.row(2) = dNd3;
    return dN;
  }

  inline Eigen::MatrixXd ShapeFunctionDerivatives(double u, double v, double w, NurbsSolid const &solid)
  {
    Eigen::MatrixXd dN(3,solid.Q.size());
    dN.row(0) = ShapeFunctionDerivative(u,v,w,0,solid);
    dN.row(1) = ShapeFunctionDerivative(u,v,w,1,solid);
    dN.row(2) = ShapeFunctionDerivative(u,v,w,2,solid);
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