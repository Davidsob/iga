#include "splines/Algorithms.h"
#include "splines/BSpline.h"
#include "splines/Nurbs.h"
#include "splines/utils/VectorOperations.h"
#include "splines/utils/Transformations.h"
#include "splines/utils/Converters.h"
#include "splines/utils/Quaternion.h"

#include "iga/ShapeFunctions.h"
#include "iga/IntegrationPoints.h"
#include "iga/SolidElementMapper.h"

#include <iostream>
#include <vector>
#include <string>

// #include <bits/stdc++.h> 
#include <cmath> 
#include <type_traits>

static std::string const python{"~/anaconda2/bin/python2.7 "};

#include <Eigen/Dense>

using namespace Eigen;

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

void sweepCurve(NurbsCurve &curve, double radius=1.0)
{
  using namespace vector_ops;
  int p = 2;
  std::vector<double> knot{0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0};
  algo::normalizeKnot(knot);
  typename NurbsCurve::matrix cpts{
    {0.0 , 0.0    , 4.0},
    {0.0 , 0.0    , 2.0},
    {0.0 , 0.0    , 0.0},
    {0.0 , 0.0    , -radius},
    {radius     , 0.0    , -radius},
  };

  // compute weights
  static double const fact{1.0/std::sqrt(2.0)};
  decltype(knot) weights(cpts.size(),1.0);
  for (size_t i = 3; i < weights.size(); i+= 2) weights[i] = fact; 

  // set properties of curve
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

void elbow(double ri, double ro, double re, NurbsSolid &solid)
{
  NurbsCurve curve; sweepCurve(curve,re);
  NurbsSurface surf; annulus(surf,ri,ro);

  solid.p = surf.p;
  solid.q = surf.q;
  solid.r = curve.p;

  solid.uknot = surf.uknot;
  solid.vknot = surf.vknot;
  solid.wknot = curve.knot;

  // now the fun part
  auto addPoints = [&solid](double cw, auto const &section)
  {
    solid.Q.insert(solid.Q.end(), section.Q.begin(), section.Q.end());
    for (auto const &w : section.weights){
      solid.weights.push_back(w*cw);
    }
  };

  {
    auto ders = spline_ops::CurveDerivatives(0,1,curve);
    auto section = transform::project(surf,{0,0,1},ders[0], ders[1]);
    addPoints(curve.weights[0],section);
  }

  {
    auto ders = spline_ops::CurveDerivatives(0.25,1,curve);
    auto section = transform::project(surf,{0,0,1},ders[0], ders[1]);
    addPoints(curve.weights[1],section);
  }

  {
    auto ders = spline_ops::CurveDerivatives(0.5,1,curve);
    auto section = transform::project(surf,{0,0,1},ders[0], ders[1]);
    addPoints(curve.weights[2],section);
  }

  {
    auto ders = spline_ops::CurveDerivatives(0.75,1,curve);
    auto section = transform::project(surf,{0,0,1},ders[0], ders[1]);
    addPoints(curve.weights[3],section);
  }

  {
    auto ders = spline_ops::CurveDerivatives(1.0,1,curve);
    auto section = transform::rotate(surf,{0,0,1}, ders[0], ders[1]);
    addPoints(curve.weights[4],section);
  }
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

void bsolid(BSplineSolid &solid)
{
  solid.p = 1;
  solid.q = 1;
  solid.r = 2;

  std::vector<double> knot{0,0,1,1};
  solid.uknot = knot;
  solid.vknot = knot;
  // solid.wknot = knot;
  solid.wknot = std::vector<double>{0,0,0,1,1,1};

  solid.Q = {
    {0,0,0},
    {1,0,0},
    {0,1,0},
    {1,1,0},

    {0,0-0.1,0.5},
    {1,0-0.1,0.5},
    {0,1-0.1,0.5},
    {1,1-0.1,0.5},

    {0,0+0.25,1},
    {1,0+0.25,1},
    {0,1+0.25,1},
    {1,1+0.25,1},
  };
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
  // std::cout << solid << std::endl;
  double u,v,w,dw;
  u = 0.8; v = 0.0; w = 0.0; dw = 1.0/20;
  std::vector<std::vector<double>> p, du, dv;
  for (int i = 0; i <= 20; i++)
  {
    p.push_back(spline_ops::SolidPoint(u,v,w,solid));
    du.push_back(spline_ops::SolidDerivative(u,v,w,1,0,solid));
    dv.push_back(spline_ops::SolidDerivative(u,v,w,1,1,solid));
    w += dw;
  }
  std::string file("output/nurbs_solid.txt");
  std::string uvec_file("output/nurbs_solid_du.txt");
  std::string vvec_file("output/nurbs_solid_dv.txt");
  spline_ops::writeVectorData(p,du, uvec_file, true, 1.0);
  spline_ops::writeVectorData(p,dv, vvec_file, true, 1.0);
  spline_ops::writeToFile(solid,file,10,1,10);
  // std::system(std::string(python + "python/plot_surface.py " + file + " ").c_str());
  std::system(std::string(python + "python/plot_surface.py " + file + " " + uvec_file + " " + vvec_file).c_str());
}

template<typename Solid>
void shapeFunctionTest(Solid const &solid)
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
  double u,v,w;
  u = 0.413; v = 0.555; w = 0.9;
  auto N = iga::ShapeFunctions(u,v,w,solid);
  std::cout << "P = " << spline_ops::SolidPoint(u,v,w,solid) << std::endl;
  auto Q = convert::to<MatrixXd>(solid.Q);
  std::cout << "P* = " << N*Q << std::endl;
  std::cout << "dP = " << spline_ops::SolidDerivative(u,v,w,1,0,solid) << std::endl;
  std::cout << "dP = " << spline_ops::SolidDerivative(u,v,w,1,1,solid) << std::endl;
  std::cout << "dP = " << spline_ops::SolidDerivative(u,v,w,1,2,solid) << std::endl;
  auto dN = iga::ShapeFunctionDerivatives(u,v,w,solid);
  std::cout << dN*Q << std::endl;
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

struct Mapper {};

void integrationPointTest()
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
  Mapper map;
  std::cout << iga::integrationPoints(map,1,1,1) << std::endl;
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

template<typename Solid>
void mapperTest(Solid const &solid)
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
  SolidElementMapper mapper(solid);
  mapper.updateElementMesh(0);

  auto ip = iga::integrationPoints(mapper,solid.p,solid.q,solid.r);
  std::cout << ip << std::endl;
  for (auto &p : ip) p.update();
  std::cout << ip << std::endl;
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

int main(int argc, char **argv)
{
  std::cout << "*** B-Spline Main ***" << std::endl;
  NurbsSolid solid; asolid(solid);
  mapperTest(solid);
  // shapeFunctionTest(solid);
  // NurbsSolid solid; elbow(0.3, 1.0, 2.0, solid);
  // integrationPointTest();
  // SolidTest(solid);
  // 
  return 0;
}
