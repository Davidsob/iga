#include "splines/Algorithms.h"
#include "splines/BSpline.h"
#include "splines/Nurbs.h"
#include "splines/SplineModifiers.h"
#include "splines/utils/VectorOperations.h"
#include "splines/utils/Transformations.h"
#include "splines/utils/Converters.h"
#include "splines/utils/Quaternion.h"

#include "iga/CurveElementMapper.h"
#include "iga/ManifoldElementMapper.h"
#include "iga/GeometricDofManager.h"
#include "iga/IntegrationPoints.h"
#include "iga/IgaIO.h"
#include "iga/ParametricMesh.h"
#include "iga/PlaneElementMapper.h"
#include "iga/Quadrature.h"
#include "iga/ShapeFunctions.h"

#include <iostream>
#include <vector>
#include <string>

// #include <bits/stdc++.h> 
#include <cmath> 
#include <type_traits>

static std::string const python{"~/anaconda2/bin/python2.7 "};

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <random>

using namespace Eigen;

void acircle(NurbsCurve &curve, double R=1.0, double th=0.0)
{
  curve.p = 2;
  curve.knot = {0,0,0,1,1,1};

  curve.Q = {
    {R,0},
    {R,R},
    {0,R},
  };

  auto const fact = std::sqrt(2)/2.0;
  curve.weights = {1.0,fact,1.0};
}

void ahemisphere(NurbsSurface &surf, double R=1.0)
{
  surf.p = 2;
  surf.q = 2;
  surf.uknot = {0,0,0,1,1,1};
  surf.vknot = {0,0,0,1,1,1};

  surf.Q = {
    {R,0,0},{R,R,0},{0,R,0},
    {R,0,0},{R,R,R},{0,R,R},
    {R,0,0},{R,0,R},{0,0,R},
  };

  auto const fact = 1.0/std::sqrt(2.0);
  surf.weights = {1.0,fact,1.0, fact,fact*fact,fact, 1.0, fact, 1.0};
}

void acircle(NurbsSurface &surf, double R)
{
  surf.p = 2;
  surf.q = 1;
  surf.uknot = {0,0,0,1,1,1};
  surf.vknot = {0,0,1,1};

  surf.Q = {
    {0,0,0},{0,0,0},{0,0,0},
    {R,0,0},{R,R,0},{0,R,0},
  };

  auto const fact = std::sqrt(2)/2.0;
  surf.weights = {1.0, fact,1.0, 1.0,fact,1.0};
}

void acylinder(NurbsSurface &surf, double R, double H)
{
  surf.p = 2;
  surf.q = 1;
  surf.uknot = {0,0,0,1,1,1};
  surf.vknot = {0,0,1,1};

  surf.Q = {
    {R,0,0},{R,R,0},{0,R,0},
    {R,0,H},{R,R,H},{0,R,H},
  };

  auto const fact = std::sqrt(2)/2.0;
  surf.weights = {1.0, fact,1.0, 1.0,fact,1.0};
}

void aplane(NurbsSurface &surf, double L, double W)
{
  surf.p = 1;
  surf.q = 1;
  surf.uknot = {0,0,1,1};
  surf.vknot = {0,0,1,1};

  surf.Q = {
    {0,0,0},{L,0,0},
    {0,W,0},{L,W,0}
  };

  surf.weights = std::vector<double>(surf.Q.size(),1);
}

void anincline(NurbsSurface &surf, double L, double W, double angle)
{
  surf.p = 1;
  surf.q = 1;
  surf.uknot = {0,0,1,1};
  surf.vknot = {0,0,1,1};

  surf.Q = {
    {0,0,0},{L*std::cos(angle),0,L*std::sin(angle)},
    {0,W,0},{L*std::cos(angle),W,L*std::sin(angle)}
  };

  surf.weights = std::vector<double>(surf.Q.size(),1);
}

void aline(NurbsCurve &curve, int p, double L=1.0, double th=0.0)
{
  curve.knot.clear();
  curve.Q.clear();

  curve.p = p;
  for (int i = 0; i <= p; i++) curve.knot.push_back(0);
  for (int i = 0; i <= p; i++) curve.knot.push_back(1);

  auto dL = L/p;
  auto C = std::cos(th);
  auto S = std::sin(th);
  for (int i = 0; i <= p; i++) curve.Q.push_back({i*dL*C,i*dL*S});
  curve.weights = std::vector<double>(curve.Q.size(),1.0);
}

void bline(NurbsCurve &curve, int p,double L = 1.0, double th=0.0)
{
  curve.knot.clear();
  curve.Q.clear();

  curve.p = p;
  for (int i = 0; i <= p; i++) curve.knot.push_back(0);
  for (int i = 0; i < p; i++)
    curve.knot.push_back(0.5);
  for (int i = 0; i <= p; i++) curve.knot.push_back(1);
  int n = curve.knot.size() - p - 1;
  auto dL = L/(n-1);
  for (int i = 0; i < n; i++) curve.Q.push_back({i*dL*std::cos(th),i*dL*std::sin(th)});
  curve.weights = std::vector<double>(curve.Q.size(),1.0);
}

void test1()
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
  double R = 1.0;
  NurbsCurve curve;
  acircle(curve,R);
  std::cout << curve << std::endl;

  size_t N = 5;
  std::vector<double> u(N+1); std::iota(u.begin(),u.end(),0); u/= N;

  auto coord = [R](double s)->std::vector<double>
  {
    return {R*std::cos(0.5*M_PI*s), R*std::sin(0.5*M_PI*s)};
  };

  auto tangent = [R](double s)->std::vector<double>
  {
    return {-0.5*M_PI*R*std::sin(0.5*M_PI*s), 0.5*M_PI*R*std::cos(0.5*M_PI*s)};
  };

  auto const Q = convert::to<Eigen::MatrixXd>(curve.Q);

  std::printf("\nPoint evaluation check....\n");
  for (auto s : u)
  {
    auto analytic = coord(s);
    auto numeric  = spline_ops::CurvePoint(s,curve);
    std::printf("\nf(%f) = ...\n",s);
    std::cout << "    analytic : " << analytic << std::endl;
    std::cout << "    numeric  : " << numeric  << std::endl;
    std::cout << "iga numeric  : " << iga::ShapeFunctions(s,curve)*Q << std::endl;

    // auto const active{iga::ActiveControlPoints(s,curve)};
    // auto const cQ = convert::to<Eigen::MatrixXd>(subVector(curve.Q,active)); // compact form
    // std::cout << "iga numeric C: " << iga::CompactShapeFunctions(s,curve)*cQ << std::endl;
  }

  std::printf("\n\nDerivative evaluation check....\n");
  for (auto s : u)
  {
    auto analytic = tangent(s);
    auto numeric  = spline_ops::CurveDerivatives(s,1,curve);
    std::printf("\nf(%f) = ...\n",s);
    std::cout << "    analytic : " << normalize(analytic) << std::endl;
    std::cout << "    numeric  : " << normalize(numeric[1] ) << std::endl;

    Eigen::MatrixXd dN = iga::ShapeFunctionDerivatives(s,curve);
    Eigen::MatrixXd iga_numeric  = dN*Q;
    std::cout << "iga numeric  : " << iga_numeric.row(0).normalized() << std::endl;

    // auto const active{iga::ActiveControlPoints(s,curve)};
    // auto const cQ = convert::to<Eigen::MatrixXd>(subVector(curve.Q,active)); // compact form
    // Eigen::MatrixXd cdN = iga::CompactShapeFunctionDerivatives(s,curve);
    // Eigen::MatrixXd ciga_numeric  = cdN*cQ;
    // std::cout << "iga numeric C: " << ciga_numeric.row(0).normalized() << std::endl;
  }

  // // plot curve for sanity
  // std::string const file = "output/arc.txt";
  // spline_ops::writeToFile(curve,file,100);
  // std::system(std::string(python + " " + "python/plot_curve.py " + file).c_str());
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

void test2()
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
  // double R = 1.0;
  // NurbsSurface surf;
  // ahemisphere(surf,R);

  // double R = 2.458e3;
  // NurbsSurface surf;
  // acircle(surf,R);
  // spline_ops::elevate(1,0,surf);
  // spline_ops::elevate(2,1,surf);

  // size_t N = 4;
  // double v = 0.876;
  // std::vector<double> u(N+1); std::iota(u.begin(),u.end(),0); u/= N;

  // auto const Q = convert::to<Eigen::MatrixXd>(surf.Q);
  // std::printf("\nPOINT EVALUATION CHECK....\n");
  // for (auto s : u)
  // {
  //   auto numeric  = spline_ops::SurfacePoint(s,v,surf);
  //   std::printf("\nf(%f,%f) = ...\n",s,v);
  //   std::cout << "    numeric  : " << numeric  << std::endl;
  //   std::cout << "iga numeric  : " << iga::ShapeFunctions(s,v,surf)*Q << std::endl;

  //   auto const active{iga::ActiveControlPoints(s,v,surf)};
  //   auto const cQ = convert::to<Eigen::MatrixXd>(subVector(surf.Q,active)); // compact form
  //   std::cout << "iga numeric C: " << iga::CompactShapeFunctions(s,v,surf)*cQ << std::endl;
  // }

  // std::printf("\n\nDERIVATIVE EVALUATION CHECK....");
  // for (auto s : u)
  // {
  //   auto numeric_ders = spline_ops::SurfaceDerivatives(s,v,2,surf); // get all derivatives up to second order
  //   Eigen::MatrixXd dN      = iga::ShapeFunctionDerivatives(s,v,surf);
  //   Eigen::RowVectorXd dNuv = iga::ShapeFunctionMixedDerivative(s,v,surf);
  //   Eigen::RowVectorXd dNuu = iga::ShapeFunctionDerivative2(s,v,0,surf);
  //   Eigen::RowVectorXd dNvv = iga::ShapeFunctionDerivative2(s,v,1,surf);
  //   Eigen::MatrixXd iga_numeric  = dN*Q;

  //   auto const active{iga::ActiveControlPoints(s,v,surf)};
  //   auto const cQ = convert::to<Eigen::MatrixXd>(subVector(surf.Q,active)); // compact form
  //   Eigen::MatrixXd cdN      = iga::CompactShapeFunctionDerivatives(s,v,surf);
  //   Eigen::RowVectorXd cdNuv = iga::CompactShapeFunctionMixedDerivative(s,v,surf);
  //   Eigen::RowVectorXd cdNuu = iga::CompactShapeFunctionDerivative2(s,v,0,surf);
  //   Eigen::RowVectorXd cdNvv = iga::CompactShapeFunctionDerivative2(s,v,1,surf);
  //   Eigen::MatrixXd ciga_numeric  = cdN*cQ;

  //   std::printf("\ndf(%f,%f) = ...\n",s,v);
  //   std::cout << "    numeric du  : " << numeric_ders[1][0] << std::endl;
  //   std::cout << "    *** iga du  : " << iga_numeric.row(0) << std::endl;
  //   std::cout << "    ***ciga du  : " << ciga_numeric.row(0) << std::endl;
  //   std::cout << "    numeric dv  : " << numeric_ders[0][1] << std::endl;
  //   std::cout << "    *** iga dv  : " << iga_numeric.row(1) << std::endl;
  //   std::cout << "    ***ciga dv  : " << ciga_numeric.row(1) << std::endl;
  //   std::cout << "    numeric duv : " << numeric_ders[1][1] << std::endl;
  //   std::cout << "    *** iga duv : " << dNuv*Q << std::endl;
  //   std::cout << "    ***ciga duv : " << cdNuv*cQ << std::endl;
  //   std::cout << "    numeric duu : " << numeric_ders[2][0] << std::endl;
  //   std::cout << "    *** iga duu : " << dNuu*Q << std::endl;
  //   std::cout << "    ***ciga duu : " << cdNuu*cQ << std::endl;
  //   std::cout << "    numeric dvv : " << numeric_ders[0][2] << std::endl;
  //   std::cout << "    *** iga dvv : " << dNvv*Q << std::endl;
  //   std::cout << "    ***ciga dvv : " << cdNvv*cQ << std::endl;
  // }

  // // // plot curve for sanity
  // std::string const file = "output/ahemisphere.txt";
  // spline_ops::writeToFile(surf,file,10,10);
  // std::system(std::string(python + " " + "python/plot_surface.py " + file).c_str());
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

template<typename Curve>
void computeArclength(Curve const &curve, double answer)
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
  // compute the length of a curve via guass quad
  // std::cout << curve << std::endl;
  // get a mapper
  CurveElementMapper mapper(curve);
  // get integration points
  auto ip = iga::integrationPoints(mapper,curve.p);
  // integrate
  auto elu = iga::meshFromSpan(curve.knot).size()-1;
  // std::printf("Curve has {%lu} elements\n",elu);
  double f = 0.0;
  double fel = 0.0;
  for (size_t i = 0; i < elu; i++)
  {
    // std::printf("\nIntegrating element(%lu)\n",i);
    // std::cout << "updating mesh..." << std::endl;
    mapper.updateElementMesh(i);
    // std::cout << "updating integration points..." << std::endl;
    for (auto &p : ip) p.update();
    // std::cout << "integrating function..." << std::endl;
    fel = 0.0;
    quadrature::gauss(ip,fel,[](auto const &p) { return 1.0; });
    f += fel;
  }
  std::printf("**** QUAD(f(x)): %10.10g\n",f);
  std::printf("**** analytic  : %10.10g\n",answer);
  std::printf("     error     : %10.10g%%\n",100.0*std::fabs(f-answer)/answer);
  GeometricDofManager::instance().clear();
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

template<typename Surface>
void computeArea(Surface const &surface, double const answer)
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
  // compute the length of a surface via guass quad
  // std::cout << surface << std::endl;
  // get a mapper
  ManifoldElementMapper mapper(surface);
  // get integration points
  auto ip = iga::integrationPoints(mapper,surface.p,surface.q);
  // std::cout << ip << std::endl;
  // integrate
  auto elu = iga::meshFromSpan(surface.uknot).size()-1;
  auto elv = iga::meshFromSpan(surface.vknot).size()-1;

  double f = 0.0;
  double fel = 0.0;
  for (size_t j = 0; j < elv; j++)
  {
    for (size_t i = 0; i < elu; i++)
    {
      // std::printf("\nIntegrating element(%lu,%lu)\n",i,j);
      // std::cout << "updating mesh..." << std::endl;
      mapper.updateElementMesh(i,j);
      std::cout << "updating integration points..." << std::endl;
      for (auto &p : ip) p.update();
      // std::cout << "integrating function..." << std::endl;
      fel = 0.0;
      quadrature::gauss(ip,fel,[](auto const &p) { return 1.0; });
      f += fel;
    }
  }
  std::printf("**** QUAD(f(x)): %10.10g\n",f);
  std::printf("**** analytic  : %10.10g\n",answer);
  std::printf("     error     : %10.10g%%\n",100.0*std::fabs(f-answer)/answer);
  GeometricDofManager::instance().clear();
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

/**
 * @brief      Computes the arclength of a curve
 */
void arclength_test()
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
  {
    std::cout << "\nCOMPUTE THE LENGTH OF A LINE:" << std::endl;
    double L = 1.13;
    double th = M_PI/3.0;
    for (int i = 1; i < 4; i++)
    {
      std::printf("Polynomial order %d:\n",i);
      NurbsCurve curve; bline(curve,i,L,th);
      computeArclength(curve,L);
    }
  }

  {
    std::cout << "\nCOMPUTE THE LENGTH OF AN ARC:" << std::endl;
    double R = 1.0;
    NurbsCurve curve;
    acircle(curve,R);
    computeArclength(curve,0.5*M_PI*R);
    for (size_t i = 0; i < 2; i++)
    {
    // spline_ops::midpointRefinement(2,curve);
      spline_ops::elevate(1,curve);
      std::printf("Polynomial order %d:\n",curve.p);
      computeArclength(curve,0.5*M_PI*R);
    }
  }
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

void surface_area_test()
{
  // {
  //   std::cout << "\nCOMPUTE THE AREA OF A PLANE:" << std::endl;
  //   double L = 10.0;
  //   double W = 0.5;
  //   NurbsSurface surf;
  //   aplane(surf,L,W);
  //   computeArea(surf,L*W);
  // }

  // {
  //   std::cout << "\nCOMPUTE THE AREA OF AN INCLINED PLANE:" << std::endl;
  //   double L = 1.13;
  //   double W = 12.15;
  //   size_t N = 5;
  //   double dth = 0.5*M_PI/N;
  //   for (size_t i = 0; i <= N; i++)
  //   {
  //     NurbsSurface surf;
  //     anincline(surf,L,W,i*dth);
  //     std::printf("Angle %f:\n",i*dth*(180.0/M_PI));
  //     computeArea(surf,L*W);
  //   }
  // }

  // {
  //   std::cout << "\nCOMPUTE THE AREA OF A CYLINDRICAL MANIFOLD:" << std::endl;
  //   double R = 2.0;
  //   double H = 0.4;
  //   NurbsSurface surf;
  //   acylinder(surf,R,H);
  //   spline_ops::elevate(1,0,surf);
  //   spline_ops::elevate(2,1,surf);
  //   double analytic = 0.5*M_PI*R*H;
  //   std::string const file = "output/surtf.txt";
  //   spline_ops::writeToFile(surf,file,30,2);
  //   std::system(std::string(python + "python/plot_surface.py " + file).c_str());
  //   // surface area test
  //   computeArea(surf,analytic);
  // }

  std::cout << "**** DEGENERATE SHAPE TESTS..." << std::endl;

  {
    std::cout << "\nCOMPUTE THE AREA OF A QUADRANT OF A CIRCLE:" << std::endl;
    double R = 3.8e3;
    NurbsSurface surf;
    acircle(surf,R);
    spline_ops::elevate(1,0,surf);
    spline_ops::elevate(2,1,surf);
    double analytic = M_PI*R*R/4.0;
    std::string const file = "output/surtf.txt";
    spline_ops::writeToFile(surf,file,20,4);
    std::system(std::string(python + "python/plot_surface.py " + file).c_str());
    // surface area test
    computeArea(surf,analytic);
  }

  {
    std::cout << "\nCOMPUTE THE AREA OF A HEMISPHERICAL MANIFOLD:" << std::endl;
    double R = 2.458;
    NurbsSurface surf;
    ahemisphere(surf,R);
    spline_ops::elevate(1,0,surf);
    spline_ops::elevate(1,1,surf);
    double analytic = 4.0*M_PI*R*R/8.0;
    std::string const file = "output/surtf.txt";
    spline_ops::writeToFile(surf,file,10,10);
    std::system(std::string(python + "python/plot_surface.py " + file).c_str());
    // surface area test
    computeArea(surf,analytic);
  }
}

int main(int argc, char **argv)
{
  // test1(); // curve derivatives
  // test2(); // surface derivatives
  // arclength_test(); // arclength test
  surface_area_test(); // surface area test
  return 0;
}
