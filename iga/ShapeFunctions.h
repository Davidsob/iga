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
  RowVectorXd ShapeFunctions(double u, double v, double w, BSplineSolid const &solid)
  {
    RowVectorXd shape(VectorXd::Zero(solid.Q.size()));

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
          shape[idx] += Nu*Nv*Nw;
        }
      }
    }
    return shape;
  }

  RowVectorXd ShapeFunctions(double u, double v, double w, NurbsSolid const &solid)
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

  std::vector<double> ShapeFunctionDerivative(double u, double v, double w, int order, int direction, BSplineSolid const &solid)
  {
    std::vector<double> dShape(solid.Q.size(),0.0);

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
              dShape[idx] += Nu*dNv*Nw;
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
              dShape[idx] += Nu*Nv*dNw;
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
              dShape[idx] += dNu*Nv*Nw;
            }
          }
        }
        break;
      }
    }

    return dShape;
  }

  std::vector<double> ShapeFunctionDerivative(double u, double v, double w, int order, int direction, NurbsSolid const &solid)
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
    // auto const shape{ShapeFunctions(u,v,w,b)};
    std::vector<double> const dshape{ShapeFunctionDerivative(u,v,w,order,direction,b)};

    // auto const wgt  = dot(shape,solid.weights);
    // auto const wgt2 = std::pow(wgt,2);
    // auto const dwgt = dot(dshape,solid.weights);

    // std::vector<double> dN{(solid.weights*(wgt*dshape - wgt*dwgt*shape)/wgt2)};
    // return dN;
    return dshape;
  }

  std::vector<std::vector<double>> ShapeFunctionDerivatives(double u, double v, double w, NurbsSolid const &solid)
  {
    std::vector<std::vector<double>> dN;
    dN.push_back(ShapeFunctionDerivative(u,v,w,1,0,solid));
    dN.push_back(ShapeFunctionDerivative(u,v,w,1,1,solid));
    dN.push_back(ShapeFunctionDerivative(u,v,w,1,2,solid));

    return dN;
  }

  Eigen::RowVectorXd parametricShapeFunction(double u, double v, double w)
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

  Eigen::MatrixXd parametricShapeFunctionDerivatives(double u, double v, double w)
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