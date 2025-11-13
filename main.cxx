
#include "splines/Algorithms.h"
#include "splines/BSpline.h"
#include "splines/Nurbs.h"
#include "splines/SplineModifiers.h"
#include "splines/DegreeReduction.h"
#include "splines/utils/VectorOperations.h"

#include "iga/IgaIO.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

#include <bits/stdc++.h> 

static std::string const python{"~/anaconda2/bin/python2.7 "};

void curveTest()
{
  using namespace vector_ops;
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  BSplineCurve curve;

  int p = 3;
  typename decltype(curve)::vector knot{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
  typename decltype(curve)::matrix cpts{
    {0.0, 0.0},
    {0.25, 1.0},
    {2.0, 0.9},
    {3.0, 0.0}
  };

  curve.p = p;
  curve.knot = knot;
  curve.Q = cpts;
  std::cout << "point on curve = " << spline_ops::CurvePoint(0.5, curve) << std::endl;
  std::cout << "derivs = " << spline_ops::CurveDerivatives(0.5, 2, curve) << std::endl;

  // std::string file("output/foo.txt");
  // spline_ops::writeToFile(curve,file);
  // std::system(std::string("~/anaconda2/bin/python2.7 python/plot_curve.py " + file).c_str());
}

void acurve(NurbsCurve &curve)
{
  int p = 3;
  typename NurbsCurve::vector knot{0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0};
  typename NurbsCurve::vector weights{1,1,1,10,20,1,1};
  typename NurbsCurve::matrix cpts{
    {0.0, 0.0},
    {0.25, 1.0},
    {0.7, 1.2},
    {0.4, -0.25},
    {1.25, -0.25},
    {1.0, 0.75},
    {1.5, 0.8}
  };

  curve.p = p;
  curve.knot = knot;
  curve.weights = weights;
  curve.Q = cpts;
}

void bcurve(BSplineCurve &curve)
{
  int p = 3;
  typename NurbsCurve::vector knot{0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0};
  typename NurbsCurve::matrix cpts{
    {0.0, 0.0},
    {0.25, 1.0},
    {0.7, 1.2},
    {0.4, -0.25},
    {1.25, -0.25},
    {1.0, 0.75},
    {1.5, 0.8}
  };

  curve.p = p;
  curve.knot = knot;
  curve.Q = cpts;
}

void asurf(BSplineSurface &surf)
{
  int p = 3; int q = 1;
  typename BSplineSurface::vector uknot{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
  typename BSplineSurface::vector vknot{0.0, 0.0, 1.0, 1.0};
  typename BSplineSurface::matrix cpts{
    {0.0, 0.0, 0.0},
    {1.0, 1.0, 0.0},
    {2.0, 1.2, 0.0},
    {3.0, 0.0, 0.0}, 

    {0.0, 0.0, 1.0},
    {1.0, 1.0, 1.0},
    {2.0, 1.2, 1.0},
    {3.0, 0.0, 1.0}
  };

  surf.p = p;
  surf.q = q;
  surf.uknot = uknot;
  surf.vknot = vknot;
  surf.Q = cpts;
}

void nurbsCurveTest()
{
  using namespace vector_ops;
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  NurbsCurve curve; acurve(curve);

  auto C = spline_ops::CurvePoint(0.5, curve);
  // make a bunch of tangents
  int N = 5;
  double du = 1.0/N; double u = 0.0;
  std::vector<std::vector<double>> pt,d;
  for (int i = 0; i <= N; i++)
  {
    auto dC = spline_ops::CurveDerivatives(u, 1, curve);
    pt.push_back(dC[0]);
    d.push_back(dC[1]);
    u += du;
  }

  std::string file("output/nurbs_curve.txt");
  std::string vec_file("output/nurbs_derivatives.txt");
  spline_ops::writeToFile(curve,file,500);
  spline_ops::writeVectorData(pt,d,vec_file,true,0.5);
  std::cout << "Plot it...." << std::endl;
  std::system(std::string(python + "python/plot_curve.py " + file + " " + vec_file).c_str());
}

void surfaceTest()
{
  using namespace vector_ops;
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  BSplineSurface surf; asurf(surf);

  std::cout << "point on surface= " << spline_ops::SurfacePoint(0.5, 0.5, surf) << std::endl;
  auto ders = spline_ops::SurfaceDerivatives(0.5, 0.5, 1, surf);
  auto Su = ders[1][0];
  auto Sv = ders[0][1];
  auto S  = ders[0][0];

  std::cout << "S = " << S << std::endl;
  std::cout << "Su = " << Su << std::endl;
  std::cout << "Sv = " << Sv << std::endl;

  std::string file("output/foo.txt");
  spline_ops::writeToFile(surf,file,8,2);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());
}

void nurbsSurfaceTest()
{
  using namespace vector_ops;
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  NurbsSurface surf;

  int p = 3;
  int q = p;
  // typename decltype(surf)::vector uknot{0.0, 0.0, 0.0, 1.0/3.0, 2.0/3.0, 1.0, 1.0, 1.0};
  typename decltype(surf)::vector uknot{0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0};
  typename decltype(surf)::vector vknot{uknot};
  typename decltype(surf)::matrix cpts{
    {1.0 , 0.0, 0.25},
    {0.6 , 0.0, 0.25},
    {0.5 , 0.0, 0.75},
    {0.25, 0.0, 0.75}, 
    {0.0 , 0.0, 1.0}, 

    {1.0 , 0.2, 0.25},
    {0.6 , 0.15, 0.25},
    {0.5 , 0.1, 0.75},
    {0.25, 0.1, 0.75}, 
    {0.0 , 0.1, 1.0}, 

    {1.0 , 0.65, -0.25},
    {0.6 , 0.6, -0.25},
    {0.5 , 0.35, 0.25},
    {0.25, 0.35, 0.25}, 
    {0.0 , 0.35, 0.5}, 

    {1.0 , 0.9, -0.25},
    {0.6 , 0.7, -0.25},
    {0.5 , 0.4, 0.25},
    {0.25, 0.4, 0.25}, 
    {0.0 , 0.4, 0.5}, 

    {1.0 , 1.3, 0.15},
    {0.6 , 1.1, 0.15},
    {0.5 , 0.6, 0.65},
    {0.25, 0.6, 0.65}, 
    {0.0 , 0.8, 0.9}
  };

  std::vector<double> weights(cpts.size(),1.0);
  auto gidx = [&uknot,p](int i, int j) {return i + j*(uknot.size()-p-1); };
  weights[gidx(1,1)] = 10.0;
  weights[gidx(1,2)] = 10.0;
  weights[gidx(2,1)] = 10.0;
  weights[gidx(2,2)] = 10.0;

  surf.p = p;
  surf.q = q;
  surf.weights = weights;
  surf.uknot = uknot;
  surf.vknot = vknot;
  surf.Q = cpts;

  auto ders = spline_ops::SurfaceDerivatives(0.7, 0.75, 1, surf);
  auto Su = ders[1][0];
  auto Sv = ders[0][1];
  auto S0  = ders[0][0];

  std::cout << "S0 = " << S0 << std::endl;
  std::cout << "Su = " << Su << std::endl;
  std::cout << "Sv = " << Sv << std::endl;

  std::string file("output/nurbs_surface.txt");
  std::string uvec_file("output/nurbs_surface_du.txt");
  std::string vvec_file("output/nurbs_surface_dv.txt");
  spline_ops::writeVectorData({S0}, {Su}, uvec_file, true, 0.2);
  spline_ops::writeVectorData({S0}, {Sv}, vvec_file, true, 0.2);
  spline_ops::writeToFile(surf,file,10,10);
  std::system(std::string(python + "python/plot_surface.py " + file + " " + uvec_file + " " + vvec_file).c_str());
}

void curveKnotInsertion()
{
  using namespace vector_ops;
  std::string file("output/nurbs_curve.txt");

  NurbsCurve curve;
  curve.p = 3;
  curve.Q = {{0,0}, {1,1}, {2,-1}, {5,-1}, {6,2}, {5,4}, {2,4},{1,1.1}};
  curve.weights = std::vector<double>(curve.Q.size(),1);
  curve.knot = {0,0,0,0,1,2,3,4,5,5,5,5};
  // acurve(curve);
  // BSplineCurve curve; bcurve(curve);
  std::cout << curve << std::endl;
  spline_ops::writeToFile(curve,file,50);
  std::system(std::string(python + "python/plot_curve.py " + file).c_str());

  std::cout << "insertion...." << std::endl;
  auto cpy = curve;
  spline_ops::midpointRefinement(3,curve);
  std::cout << curve << std::endl;
  spline_ops::writeToFile(curve,file,50);
  std::system(std::string(python + "python/plot_curve.py " + file).c_str());
}

void curveElevation()
{
  using namespace vector_ops;
  std::string file("output/nurbs_curve.txt");

  NurbsCurve curve;
  curve.p = 3;
  curve.Q = {{0,0}, {1,1}, {2,-1}, {5,-1}, {6,2}, {5,4}, {2,4},{1,1.1}};
  curve.weights = std::vector<double>(curve.Q.size(),1);
  curve.knot = {0,0,0,0,1,2,3,4,5,5,5,5};
  // acurve(curve);
  // BSplineCurve curve; bcurve(curve);
  std::cout << curve << std::endl;
  // spline_ops::writeToFile(curve,file,50);
  // std::system(std::string(python + "python/plot_curve.py " + file).c_str());

  std::cout << "elevate...." << std::endl;
  auto cpy = curve;
  spline_ops::elevate(1,cpy);
  std::cout << cpy << std::endl;


  // compare point by point
  size_t N = 100;
  std::vector<double> iso(N+1); std::iota(iso.begin(),iso.end(),0); iso/= N;
  std::vector<double> isoline;
  std::vector<double> error;
  for (auto x : iso)
  {
    auto p = spline_ops::CurvePoint(x,curve);
    auto pe = spline_ops::CurvePoint(x,cpy);
    error.push_back(norm(p-pe));
  }

  std::string isofile("output/elevation_error.txt");
  IO::writeXYdata(iso,error,isofile);
  std::system(std::string(python + "python/plot_xy.py " + isofile).c_str());
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;

  // see elevated curve
  spline_ops::writeToFile(cpy,file,50);
  std::system(std::string(python + "python/plot_curve.py " + file).c_str());
}

void aline(BSplineCurve &curve, int p)
{
  curve.p = p;
  for (int i = 0; i <= p; i++) curve.knot.push_back(0);
  for (int i = 0; i <= p; i++) curve.knot.push_back(1);
  for (int i = 0; i <= p; i++) curve.Q.push_back({double(i),0});
}

void aline2(BSplineCurve &curve, int p)
{
  curve.p = p;
  for (int i = 0; i <= p; i++) curve.knot.push_back(0);
  // for (int i = 0; i < p; i++)
  //   curve.knot.push_back(0.5);
  for (int i = 0; i <= p; i++) curve.knot.push_back(1);
  int n = curve.knot.size() - p - 1;
  for (int i = 0; i < n; i++) curve.Q.push_back(std::vector<double>({double(i),0.0}));
}


void curveReduction()
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
  using namespace vector_ops;
  std::string file("output/nurbs_curve.txt");

  // NurbsCurve curve;
  // curve.p = 3;
  // curve.Q = {{0,0}, {1,1}, {2,-1}, {5,-1}, {6,2}, {5,4}, {2,4},{1,1.1}};
  // curve.weights = std::vector<double>(curve.Q.size(),1);
  // curve.knot = {0,0,0,0,1,2,3,4,5,5,5,5};
  // 
  BSplineCurve curve; aline2(curve,10);
  std::cout << curve << std::endl;
  spline_ops::writeToFile(curve,file,50);
  std::system(std::string(python + "python/plot_curve.py " + file).c_str());

  // reduce
  std::cout << "reduce ...." << std::endl;
  spline_ops::reduce(curve,1e-4);
  std::cout << curve << std::endl;
  spline_ops::writeToFile(curve,file,50);
  std::system(std::string(python + "python/plot_curve.py " + file).c_str());
  // 
  // compare point by point
  // size_t N = 100;
  // std::vector<double> iso(N+1); std::iota(iso.begin(),iso.end(),0); iso/= N;
  // std::vector<double> isoline;
  // std::vector<double> error;
  // for (auto x : iso)
  // {
  //   auto p = spline_ops::CurvePoint(x,curve);
  //   auto pe = spline_ops::CurvePoint(x,cpy);
  //   error.push_back(norm(p-pe));
  // }

  // std::string isofile("output/elevation_error.txt");
  // IO::writeXYdata(iso,error,isofile);
  // std::system(std::string(python + "python/plot_xy.py " + isofile).c_str());
  // std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

void surfaceKnotInsertion()
{
  using namespace vector_ops;
  std::string file("output/nurbs_surface.txt");
  BSplineSurface surf; asurf(surf);
  // BSplineCurve curve; bcurve(curve);
  std::cout << surf << std::endl;
  spline_ops::writeToFile(surf,file,5,5);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());

  std::cout << "insertion...." << std::endl;
  spline_ops::insertKnot(0.3,1,0,surf);
  std::cout << surf << std::endl;
  spline_ops::writeToFile(surf,file,20,5);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());
}

void surfaceElevation()
{
  using namespace vector_ops;
  std::string file("output/nurbs_surface.txt");
  BSplineSurface surf; asurf(surf);
  // BSplineCurve curve; bcurve(curve);
  spline_ops::writeToFile(surf,file,5,5);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());

  std::cout << "insertion...." << std::endl;
  BSplineSurface cpy = surf;
  spline_ops::elevate(1,1,cpy);
  std::cout << surf << std::endl;
  std::cout << cpy << std::endl;
  spline_ops::writeToFile(cpy,file,20,5);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());

  size_t N = 20;
  std::vector<double> u(N+1); std::iota(u.begin(),u.end(),0); u/= N;
  std::vector<double> v(N+1); std::iota(v.begin(),v.end(),0); v/= N;
  std::vector<double> error;
  for (auto x : u)
  {
    for (auto y : v)
    {
      auto p  = spline_ops::SurfacePoint(x,y,surf);
      auto pe = spline_ops::SurfacePoint(x,y,cpy);
      error.push_back(norm(p-pe));
    }
  }

  std::cout << "error = " << norm(error)/error.size() << std::endl;
}

void surfaceReduction()
{
  using namespace vector_ops;
  std::string file("output/nurbs_surface.txt");
  BSplineSurface surf; asurf(surf);
  std::cout << surf << std::endl;
  spline_ops::writeToFile(surf,file,5,5);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());

  std::cout << "insertion...." << std::endl;
  BSplineSurface cpy = surf;
  spline_ops::elevate(1,1,cpy);
  std::cout << cpy << std::endl;
  spline_ops::writeToFile(cpy,file,20,5);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());

  // reduce
  spline_ops::reduce(0,cpy,1);
  std::cout << cpy << std::endl;
  spline_ops::writeToFile(cpy,file,20,5);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());

  size_t N = 20;
  std::vector<double> u(N+1); std::iota(u.begin(),u.end(),0); u/= N;
  std::vector<double> v(N+1); std::iota(v.begin(),v.end(),0); v/= N;
  std::vector<double> error;
  for (auto x : u)
  {
    for (auto y : v)
    {
      auto p  = spline_ops::SurfacePoint(x,y,surf);
      auto pe = spline_ops::SurfacePoint(x,y,cpy);
      error.push_back(norm(p-pe));
    }
  }

  std::cout << "error = " << norm(error)/error.size() << std::endl;
}

int main(int argc, char **argv)
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
  // curveReduction();
  surfaceReduction();
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
  return 0;
}
