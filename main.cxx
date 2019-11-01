
#include "splines/Algorithms.h"
#include "splines/BSpline.h"
#include "splines/Nurbs.h"
#include "splines/utils/VectorOperations.h"

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

void nurbsCurveTest()
{
  using namespace vector_ops;
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  NurbsCurve curve;

  int p = 3;
  typename decltype(curve)::vector knot{0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0};
  typename decltype(curve)::vector weights{1,1,1,10,20,1,1};
  typename decltype(curve)::matrix cpts{
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
  BSplineSurface surf;

  int p = 3; int q = 1;
  typename decltype(surf)::vector uknot{0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
  typename decltype(surf)::vector vknot{0.0, 0.0, 1.0, 1.0};
  typename decltype(surf)::matrix cpts{
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
  std::cout << "### (" << __LINE__ << ")" << std::endl;
  auto Sw = spline_ops::SurfacePoint(0.5, 0.5, surf);
  std::cout << "### (" << __LINE__ << ")" << std::endl;
  std::cout << "Point on surface = " << Sw << std::endl;
  // auto ders = spline_ops::SurfaceDerivatives(0.5, 0.5, 1, surf);
  // auto Su = ders[1][0];
  // auto Sv = ders[0][1];
  // auto S  = ders[0][0];

  // std::cout << "S = " << S << std::endl;
  // std::cout << "Su = " << Su << std::endl;
  // std::cout << "Sv = " << Sv << std::endl;

  std::string file("output/nurbs_surface.txt");
  spline_ops::writeToFile(surf,file,50,50);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());
}

int main(int argc, char **argv)
{
  std::cout << "*** B-Spline Main ***" << std::endl;
  // nurbsCurveTest();
  nurbsSurfaceTest();
  return 0;
}
