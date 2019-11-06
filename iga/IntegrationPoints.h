#pragma once

#include <cassert>
#include <iostream>

struct IntegrationPoint
{
  int id;
  double r,s,t,jac,w;
  friend std::ostream &operator<<(std::ostream &os, IntegrationPoint const &ip);
};

std::ostream &operator<<(std::ostream &os, IntegrationPoint const &ip)
{
  os << "Integration point(" << ip.id << "):" << std::endl;
  os << "*** weight   = " << ip.w << std::endl;
  os << "*** jacobian = " << ip.jac << std::endl;
  os << "*** gp       = {" << ip.r << ", " << ip.s << ", " << ip.t << "}" << std::endl;
  return os;
}

namespace iga
{
  std::vector<double> gaussLegendrePoints(int order);
  std::vector<double> gaussLegendreWeights(int order);

  std::vector<IntegrationPoint> _integrationPoints(int order)
  {
    auto const gp{gaussLegendrePoints(order)};
    auto const w{gaussLegendreWeights(order)};
    std::vector<IntegrationPoint> points(gp.size(),IntegrationPoint());
    for (size_t i = 0; i < gp.size(); i++)
    {
      points[i].id = i;
      points[i].r = gp[i];
      points[i].w = w[i];
    }
    return points;
  }

  std::vector<IntegrationPoint> _integrationPoints(int order1, int order2)
  {
    auto const gp1{gaussLegendrePoints(order1)};
    auto const gp2{gaussLegendrePoints(order2)};
    auto const w1{gaussLegendreWeights(order1)};
    auto const w2{gaussLegendreWeights(order2)};
    std::vector<IntegrationPoint> points(gp1.size()*gp2.size(),IntegrationPoint());
    size_t k = 0;
    for (size_t j = 0; j < gp2.size(); j++)
    {
      for (size_t i = 0; i < gp1.size(); i++)
      {
        points[k].id = k;
        points[k].r = gp1[i];
        points[k].s = gp2[j];
        points[k].w = w1[i]*w2[j];
        k++;
      }
    }
    return points;
  }

  std::vector<IntegrationPoint> _integrationPoints(int order1, int order2, int order3)
  {
    auto const gp1{gaussLegendrePoints(order1)};
    auto const gp2{gaussLegendrePoints(order2)};
    auto const gp3{gaussLegendrePoints(order3)};
    auto const w1{gaussLegendreWeights(order1)};
    auto const w2{gaussLegendreWeights(order2)};
    auto const w3{gaussLegendreWeights(order3)};
    std::vector<IntegrationPoint> points(gp1.size()*gp2.size()*gp3.size(),IntegrationPoint());
    size_t kk = 0;
    for (size_t k = 0; k < gp3.size(); k++)
    {
      for (size_t j = 0; j < gp2.size(); j++)
      {
        for (size_t i = 0; i < gp1.size(); i++)
        {
          points[kk].id = kk;
          points[kk].r = gp1[i];
          points[kk].s = gp2[j];
          points[kk].t = gp3[k];
          points[kk].w = w1[i]*w2[j]*w3[k];
          kk++;
        }
      }
    }
    return points;
  }

  template<typename ...Args>
  std::vector<IntegrationPoint> integrationPoints(Args const ...args)
  {
    return _integrationPoints(args...);
  }

  std::vector<double> gaussLegendrePoints(int order)
  {
    std::vector<double> pts;
    switch(order)
    {
      case 0:
      {
        pts.push_back(0);
      } break;

      case 1:
      {
        static const double x{1.0/std::sqrt(3.0)};
        pts.push_back(-x);
        pts.push_back( x);
      } break;

      case 2:
      {
        static const double x{std::sqrt(3.0/5.0)};
        pts.push_back(-x);
        pts.push_back(0.0);
        pts.push_back( x);
      } break;

      case 3:
      {
        static const double x1{std::sqrt(3.0/7.0 - (2.0/7.0)*std::sqrt(6.0/5.0))};
        static const double x2{std::sqrt(3.0/7.0 + (2.0/7.0)*std::sqrt(6.0/5.0))};
        pts.push_back(-x2);
        pts.push_back(-x1);
        pts.push_back( x1);
        pts.push_back( x2);
      } break;

      case 4:
      {
        static const double x1{std::sqrt(5.0 - 2.0*std::sqrt(10.0/7.0))};
        static const double x2{std::sqrt(5.0 + 2.0*std::sqrt(10.0/7.0))};
        pts.push_back(-x2);
        pts.push_back(-x1);
        pts.push_back(0.0);
        pts.push_back( x1);
        pts.push_back( x2);
      } break;

      default: assert(false); break;
    }

    return pts;
  }

  std::vector<double> gaussLegendreWeights(int order)
  {
    std::vector<double> pts;
    switch(order)
    {
      case 0:
      {
        pts.push_back(2.0);
      } break;

      case 1:
      {
        pts.push_back(1.0);
        pts.push_back(1.0);
      } break;

      case 2:
      {
        static const double x1{5.0/9.0};
        static const double x2{8.0/9.0};
        pts.push_back(x1);
        pts.push_back(x2);
        pts.push_back(x1);
      } break;

      case 3:
      {
        static const double x1{(18.0 + std::sqrt(30.0))/36.0};
        static const double x2{(18.0 - std::sqrt(30.0))/36.0};
        pts.push_back(x2);
        pts.push_back(x1);
        pts.push_back(x1);
        pts.push_back(x2);
      } break;

      case 4:
      {
        static const double x1{(322.0 + 13.0*std::sqrt(70.0))/900.0};
        static const double x2{(322.0 - 13.0*std::sqrt(70.0))/900.0};
        static const double x3{128.0/255.0};
        pts.push_back(x2);
        pts.push_back(x1);
        pts.push_back(x3);
        pts.push_back(x1);
        pts.push_back(x2);
      } break;

      default: assert(false); break;
    }

    return pts;
  }
}