#include "splines/Algorithms.h"
#include "splines/BSpline.h"
#include "splines/Nurbs.h"
#include "splines/utils/VectorOperations.h"

#include <iostream>
#include <vector>
#include <string>

#include <bits/stdc++.h> 
#include <type_traits>

static std::string const python{"~/anaconda2/bin/python2.7 "};

void circle(NurbsCurve &curve, double radius=1.0)
{
  using namespace vector_ops;
  int p = 2;
  std::vector<double> knot{0.0, 0.0, 0.0, 1, 1, 2, 2, 3, 3, 4, 4, 4};
  algo::normalizeKnot(knot);
  typename NurbsCurve::matrix cpts{
    {radius , 0.0    , 0.0},
    {radius , radius , 0.0},
    {0.0    , radius , 0.0},
    {-radius, radius , 0.0},
    {-radius, 0.0    , 0.0},
    {-radius, -radius, 0.0},
    {0.0    , -radius, 0.0},
    {radius , -radius, 0.0},
    {radius , 0.0    , 0.0},
  };

  static double const fact{1.0/std::sqrt(2.0)};
  decltype(knot) weights(cpts.size(),1);
  for (size_t i = 1; i < weights.size(); i+= 2)
    weights[i] = fact; 
  curve.p = p;
  curve.knot = knot;
  curve.weights = weights;
  curve.Q = cpts;
}

void annulus(NurbsSurface &surf, double ri, double ro)
{
  using namespace spline_ops;

  NurbsCurve ci, co;
  circle(ci,ri);
  circle(co,ro);

  surf.p = 2;
  surf.q = 1;
  surf.Q = co.Q;
  surf.Q.insert(surf.Q.end(), ci.Q.begin(), ci.Q.end());
  surf.weights = ci.weights;
  surf.weights.insert(surf.weights.end(), co.weights.begin(), co.weights.end());

  surf.uknot = ci.knot;
  surf.vknot = std::vector<double>({0,0,1,1});
}

void asolid(NurbsSolid &solid)
{
  solid.p = 1;
  solid.q = 1;
  solid.r = 2;

  std::vector<double> knot{0,0,1,1};
  solid.uknot = knot;
  solid.vknot = knot;
  solid.wknot = std::vector<double>{0,0,0,1,1,1};

  solid.Q = {
    {0,0,0},
    {1,0,0},
    {0,1,0},
    {1,1,0},

    {0,0-0.2,0.95},
    {1,0-0.2,0.95},
    {0,1-0.2,0.95},
    {1,1-0.2,0.95},

    {0,0+0.2,1},
    {1,0+0.2,1},
    {0,1+0.2,1},
    {1,1+0.2,1},
  };

  solid.weights = std::vector<double>(solid.Q.size(),1);
}

template<typename Curve>
void CurveTest(Curve const &curve)
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  using namespace vector_ops;
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
  spline_ops::writeToFile(curve,file,50);
  spline_ops::writeVectorData(pt,d,vec_file,true,0.5);
  std::cout << "Plot it...." << std::endl;
  std::system(std::string(python + "python/plot_curve.py " + file + " " + vec_file).c_str());
}


template<typename Surface>
void SurfaceTest(Surface const &surf)
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  using namespace vector_ops;
  std::cout << spline_ops::SurfacePoint(0.5, 0.5, surf) << std::endl;
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
  spline_ops::writeToFile(surf,file,50,1);
  std::system(std::string(python + "python/plot_surface.py " + file + " " + uvec_file + " " + vvec_file).c_str());
}

template<typename Solid>
void SolidTest(Solid const &solid)
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  using namespace vector_ops;
  std::cout << solid << std::endl;
  std::cout << spline_ops::SolidPoint(0.1, 0.9, 0.5, solid) << std::endl;
  std::cout << spline_ops::SolidPoint(0.5, 0.5, 0.0, solid) << std::endl;
  std::cout << spline_ops::SolidPoint(0.78, 1.0, 1.0, solid) << std::endl;

  std::string file("output/nurbs_solid.txt");
  spline_ops::writeToFile(solid,file,2,2,5);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());
}

int main(int argc, char **argv)
{
  std::cout << "*** B-Spline Main ***" << std::endl;
  // NurbsCurve curve; circle(curve);
  // CurveTest(curve);
  // NurbsSurface surface; annulus(surface,1.0,1.25);
  // SurfaceTest(surface);
  NurbsSolid solid; asolid(solid);
  SolidTest(solid);
  return 0;
}
