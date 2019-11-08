#include "splines/Algorithms.h"
#include "splines/BSpline.h"
#include "splines/Nurbs.h"
#include "splines/utils/VectorOperations.h"
#include "splines/utils/Transformations.h"
#include "splines/utils/Converters.h"
#include "splines/utils/Quaternion.h"

#include "iga/ShapeFunctions.h"
#include "iga/IntegrationPoints.h"
#include "iga/IgaIO.h"
#include "iga/ParametricMesh.h"
#include "iga/SolidElementMapper.h"
#include "iga/Quadrature.h"

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

void arc(NurbsCurve &curve, double radius=1.0)
{
  using namespace vector_ops;
  int p = 2;
  std::vector<double> knot{0.0, 0.0, 0.0, 1, 1, 2,2, 3,3,3};
  algo::normalizeKnot(knot);
  typename NurbsCurve::matrix cpts{
    {radius , 0.0    , 0.0},
    {radius , radius , 0.0},
    {0.0    , radius , 0.0},
    {-radius, radius , 0.0},
    {-radius, 0.0    , 0.0},
    {-radius, -radius, 0.0},
    {0.0, -radius, 0.0},
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

void sector(NurbsSurface &surf, double ri, double ro)
{
  using namespace spline_ops;

  NurbsCurve ci, co;
  arc(ci,ri);
  arc(co,ro);

  surf.p = 2;
  surf.q = 1;
  surf.Q = ci.Q;
  surf.Q.insert(surf.Q.end(), co.Q.begin(), co.Q.end());
  surf.weights = ci.weights;
  surf.weights.insert(surf.weights.end(), co.weights.begin(), co.weights.end());

  surf.uknot = co.knot;
  surf.vknot = std::vector<double>({0,0,1,1});
}

void annulus(NurbsSurface &surf, double ri, double ro)
{
  using namespace spline_ops;

  NurbsCurve ci, co;
  circle(ci,ri);
  circle(co,ro);

  surf.p = 2;
  surf.q = 1;
  surf.Q = ci.Q;
  surf.Q.insert(surf.Q.end(), co.Q.begin(), co.Q.end());
  surf.weights = ci.weights;
  surf.weights.insert(surf.weights.end(), co.weights.begin(), co.weights.end());

  surf.uknot = co.knot;
  surf.vknot = std::vector<double>({0,0,1,1});
}

void pipe(double ri, double ro, NurbsSolid &solid)
{
  NurbsSurface surf; annulus(surf,ri,ro);

  solid.p = surf.p;
  solid.q = surf.q;
  solid.r = 2;

  solid.uknot = surf.uknot;
  solid.vknot = surf.vknot;
  solid.wknot = {0,0,0,1,1,1};

  // now the fun part
  auto addPoints = [&solid](double cw, auto const &section)
  {
    solid.Q.insert(solid.Q.end(), section.Q.begin(), section.Q.end());
    for (auto const &w : section.weights){
      solid.weights.push_back(w*cw);
    }
  };

  addPoints(1.0, surf);// bottom of pipe

  {
    auto section = transform::translate(surf,{0,0,2});
    addPoints(1,section);
  }

  {
    auto section = transform::translate(surf,{0,0,4});
    addPoints(1,section);
  }
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
  solid.wknot = std::vector<double>{0,0,0,0.5,0.5,1,1,1};

  solid.Q = {
    {0,0,0},
    {1,0,0},
    {0,1,0},
    {1,1,0},

    {0,0,0.9},
    {1,0,0.9},
    {0,1,0.9},
    {1,1,0.9},

    {0,0,1},
    {1,0,1},
    {0,1,1},
    {1,1,1},

    {0,0,1.4},
    {1,0,1.4},
    {0,1,1.4},
    {1,1,1.4},

    {0,0,2},
    {1,0,2},
    {0,1,2},
    {1,1,2},
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
  std::string file("output/nurbs_solid.txt");
  spline_ops::writeToFile(solid,file,10,1,10);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());
}

struct Plane
{
  std::vector<double> x;
  std::vector<double> n;
};

template<typename Shape>
std::vector<size_t> pointsOnPlane(Shape const &shape, Plane const &plane)
{
  auto on_plane = [&plane](auto const &p)
  {
    return algo::equal(dot(plane.n,p-plane.x),0.0);
  };

  std::vector<size_t> found;
  for (size_t i = 0; i < shape.Q.size(); i++)
  {
    if (on_plane(shape.Q[i])) found.push_back(i);
  }

  return found;
}

template<typename Solid>
void heatCondution(Solid const &solid)
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
// #################################
// Weak forms 
// #################################
  auto stiffness = [](auto const &p)->Eigen::MatrixXd
  {
    auto const grad = p.mapper.grad(p.para);
    auto kel = grad.transpose()*grad;
    return kel;
  };

  // Compute stiffness
  size_t dof = solid.Q.size();
  Eigen::MatrixXd K(Eigen::MatrixXd::Zero(dof,dof));
  SolidElementMapper mapper(solid);
  auto ip = iga::integrationPoints(mapper,solid.p,solid.q,solid.r);

  auto elu = iga::meshFromSpan(solid.uknot).size()-1;
  auto elv = iga::meshFromSpan(solid.vknot).size()-1;
  auto elw = iga::meshFromSpan(solid.wknot).size()-1;
  std::printf("Default mesh has {%lu, %lu, %lu} elements\n",elu,elv,elw);
  for (size_t k = 0; k < elw; k++)
    for (size_t j = 0; j < elv; j++)
      for (size_t i = 0; i < elu; i++)
      {
        std::printf("\nIntegrating element(%lu,%lu,%lu)\n",i,j,k);
        mapper.updateElementMesh(i,j,k);
        std::cout << "updating integration points..." << std::endl;
        for (auto &p : ip) p.update();
        quadrature::gauss(ip,K,stiffness);
      }

// #################################
// Boundary conditions
// #################################
  // compute constraints
  double Thot,Tcold;
  Thot = 500; Tcold = 200;
  Plane plane;
  plane.x = {0,0,4};
  plane.n = {0,0,1};
  // x- plane
  auto bc_hot = pointsOnPlane(solid,plane);
  // x+ plane
  // plane.x = {3,0,0};
  // plane.n = {1,0,0};
  plane.x = {0,0,0};
  plane.n = {0,0,-1};
  auto bc_cold = pointsOnPlane(solid,plane);

  std::printf("Applying boundary conditions for hot surface to %lu cpts\n",bc_hot.size());
  std::printf("Applying boundary conditions for cold surface to %lu cpts\n",bc_cold.size());

  auto C = bc_hot.size() + bc_cold.size();
  Eigen::MatrixXd Kc(C,dof); Kc.setZero();
  Eigen::VectorXd Rc(C); Rc.setZero();

  size_t k = 0;
  for (auto idx : bc_hot)
  {
    Kc(k,idx) = 1;
    Rc[k] = Thot;
    k++;
  }

  for (auto idx : bc_cold)
  {
    Kc(k,idx) = 1;
    Rc[k] = Tcold;
    k++;
  }

// #################################
// compose matrix 
// #################################
  Eigen::MatrixXd A(dof+C,dof+C); A.setZero();
  Eigen::VectorXd F(dof+C); F.setZero();

  F.tail(C) = Rc;
  A.block(0,0,dof,dof) = K;
  A.block(dof,0,C,dof) = Kc;
  A.block(0,dof,dof,C) = Kc.transpose();
// #################################
// solve 
// #################################
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A.sparseView());
  if(solver.info()!=Eigen::ComputationInfo::Success)
  {
    std::cout << "decomposition failed" << std::endl;
    return;
  }

  Eigen::VectorXd const T = solver.solve(F).head(dof);
  if(solver.info()!=Eigen::ComputationInfo::Success)
  {
    std::cout << "Solve failed" << std::endl;
    return;
  }


// #################################
// write solution 
// #################################
  // std::string file("output/nurbs_heat_conduction.txt");
  // IO::writeSolutionToFile(solid,convert::to<std::vector<double>>(T),file,30,5,50);
  // std::system(std::string(python + "python/plot_surface.py " + file).c_str());

  // check along u=v=const
  NurbsSolid wsolution; IO::geometryWithSolution(solid,convert::to<std::vector<double>>(T),wsolution);
  size_t N = 100;
  std::vector<double> w(N+1); std::iota(w.begin(),w.end(),0); w/= N;
  double u = 0.5; double v = 1.0;
  std::vector<double> isoline;
  for (auto x : w)
  {
    auto s = spline_ops::SolidPoint(u,v,x,wsolution);
    isoline.push_back(s.back());
  }

  std::string isofile("output/nurbs_heat_iso.txt");
  IO::writeXYdata(w,isoline,isofile);
  std::system(std::string(python + "python/plot_xy.py " + isofile).c_str());
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

template<typename Solid>
void shapeFunctionTest(Solid const &solid, size_t n)
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1

  auto rand_point = [&dis,&gen]()->std::vector<double>
  { return {dis(gen),dis(gen),dis(gen)}; };

  auto const Q = convert::to<Eigen::MatrixXd>(solid.Q);
  for (size_t i = 0; i < n; i++) {
    auto p = rand_point();
    auto S1 = spline_ops::SolidPoint(p[0],p[1],p[2],solid);
    auto N = iga::ShapeFunctions(p[0],p[1],p[2],solid);
    auto S2 = N*Q;
    auto dS = S2-convert::to<RowVectorXd>(S1);

    if (!algo::equal(dS.norm(),0.0))
      std::cout << "POINT EVALUATION FAILED!" << std::endl;
    if (!algo::equal(N.sum(),1.0))
      std::cout << "PARTITION OF UNITY FAILED!" << std::endl;
  }
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

template<typename Solid>
void derivativeTest(Solid const &solid, size_t n)
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1

  auto rand_point = [&dis,&gen]()->std::vector<double>
  { return {dis(gen),dis(gen),dis(gen)}; };

  auto const Q = convert::to<Eigen::MatrixXd>(solid.Q);
  for (size_t i = 0; i < n; i++) {
    auto p = rand_point();
    auto Su = spline_ops::SolidDerivative(p[0],p[1],p[2],1,0,solid);
    auto Sv = spline_ops::SolidDerivative(p[0],p[1],p[2],1,1,solid);
    auto Sw = spline_ops::SolidDerivative(p[0],p[1],p[2],1,2,solid);
    std::vector<std::vector<double>> tmp{Su,Sv,Sw};
    auto grad1 = convert::to<Eigen::MatrixXd>(tmp);

    auto dN = iga::ShapeFunctionDerivatives(p[0],p[1],p[2],solid);
    auto grad = dN*Q;

    auto dg = grad-grad1;
    if (!algo::equal(dg.norm(),0.0))
      std::cout << "POINT EVALUATION FAILED!" << std::endl;
  }
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

int main(int argc, char **argv)
{
  NurbsSolid solid;
  // bsolid(solid);
  pipe(1,2,solid);
  // elbow(1,2,3,solid);
  // shapeFunctionTest(solid,10000);
  // derivativeTest(solid,25);
  heatCondution(solid);
  // SolidTest(solid);
  return 0;
}
