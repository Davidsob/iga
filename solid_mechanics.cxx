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

static const Eigen::IOFormat CleanFormat(4, 0, ", ", "\n", "[", "]");

void arc(NurbsCurve &curve, double radius=1.0)
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
    {radius*std::cos(0.01) , -radius*std::sin(0.01)    , 0.0},
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
  NurbsSurface surf; sector(surf,ri,ro);

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
void simulation(Solid const &solid)
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
// #################################
// Operators 
// #################################
  auto Bop = [](auto const &p)->Eigen::MatrixXd
  {
    auto const grad = p.mapper.grad(p.para);
    auto const dof = grad.cols();
    Eigen::MatrixXd b(6,3*dof); b.setZero();

    for (int i = 0; i < dof; i++)
    {
      b(0,i*3)   = grad(0,i);
      b(1,i*3+1) = grad(1,i);
      b(2,i*3+2) = grad(2,i);

      b(3,i*3)   = grad(1,i);
      b(3,i*3+1) = grad(0,i);

      b(4,i*3+1) = grad(2,i);
      b(4,i*3+2) = grad(1,i);

      b(5,i*3)   = grad(2,i);
      b(5,i*3+2) = grad(1,i);
    }
    return b;
  };

  auto Nop = [](auto const &p)->Eigen::MatrixXd
  {
    auto const shape = p.mapper.shape(p.para);
    auto const dof = shape.size();
    Eigen::MatrixXd N(3,3*dof); N.setZero();

    for (int i = 0; i < dof; i++)
    {
      N(0,i*3)   = shape[i];
      N(1,i*3+1) = shape[i];
      N(2,i*3+2) = shape[i];
    }
    return N;
  };

// #################################
// Material 
// #################################
  double const E = 70e9;
  double const nu = 0.29;
  double const rho = 2700.0;
  auto matl = [E,nu](auto const &p)->Eigen::MatrixXd
  {
    Eigen::MatrixXd D(6,6); D.setZero();
    double c = (1.0 - 2.0*nu);
    double a = E/((1.0 + nu)*c);
    double b = a*(1.0 - nu);
    double d = a*nu;
    c *= (0.5 *a);

    D(0,0) = b;
    D(0,1) = d;
    D(0,2) = d;

    D(1,0) = d;
    D(1,1) = b;
    D(1,2) = d;

    D(2,0) = d;
    D(2,1) = d;
    D(2,2) = b;

    D(3,3) = c;
    D(4,4) = c;
    D(5,5) = c;

    return D;
  };

// #################################
// Weak forms 
// #################################
  Vector3d accel; accel << 0.0, 0.0, -9.81;
  auto gravity = [&Nop, &rho, &accel](auto const &p)->Eigen::VectorXd
  {
    auto const n = Nop(p);
    return rho*(n.transpose()*accel);
  };

  auto stiffness = [&Bop,&matl](auto const &p)->Eigen::MatrixXd
  {
    auto const b = Bop(p);
    auto const d = matl(p);
    return b.transpose()*(d*b);
  };

  // Compute stiffness
  size_t dof = 3*solid.Q.size();
  Eigen::MatrixXd K(dof,dof); K.setZero();
  Eigen::VectorXd f(dof); f.setZero();
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
        quadrature::gauss(ip,f,gravity);
      }

// #################################
// Boundary conditions
// #################################
  // compute constraints
  // fix at the top
  Plane plane;
  plane.x = {0,0,4};
  plane.n = {0,0,1};
  // x- plane
  auto fixed = pointsOnPlane(solid,plane);

  std::printf("Applying boundary conditions for fixed surface to %lu cpts\n",fixed.size());

  auto C = 3*fixed.size();
  Eigen::MatrixXd Kc(C,dof); Kc.setZero();
  Eigen::VectorXd Rc(C); Rc.setZero();

  size_t k = 0;
  for (auto const &idx : fixed)
  {
    Kc(k++,3*idx)   = 1;
    Kc(k++,3*idx+1) = 1;
    Kc(k++,3*idx+2) = 1;
  }

// #################################
// compose matrix 
// #################################
  Eigen::MatrixXd A(dof+C,dof+C); A.setZero();
  Eigen::VectorXd F(dof+C); F.setZero();
  F << f,Rc;
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

  Eigen::VectorXd const u = solver.solve(F).head(dof);
  if(solver.info()!=Eigen::ComputationInfo::Success)
  {
    std::cout << "Solve failed" << std::endl;
    return;
  }

// #################################
// write solution 
// #################################
  // reshape solution
  auto const uvw = Eigen::Map<const Eigen::MatrixXd>(u.data(),3,solid.Q.size()).transpose();
  // compute disp mag and deform the mesh
  double scale = 1e5;
  Solid deformed = solid;
  Eigen::VectorXd disp(solid.Q.size());
  for (size_t i = 0; i < solid.Q.size(); i++)
  {
    auto const x = convert::to<std::vector<double>>(Eigen::RowVectorXd(uvw.row(i)));
    disp[i] = norm(x);
    deformed.Q[i] += scale*x;
  }

  std::string file("output/nurbs_stress.txt");
  IO::writeSolutionToFile(deformed,uvw.col(2),file,20,10,20);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());

  // check along v=w=const
  NurbsSolid wsolution; IO::geometryWithSolution(solid,uvw.col(2),wsolution);
  size_t N = 100;
  std::vector<double> iso(N+1); std::iota(iso.begin(),iso.end(),0); iso/= N;
  double v = 1.0; double w = 1.0;
  std::vector<double> isoline;
  for (auto x : iso)
  {
    auto s = spline_ops::SolidPoint(x,v,w,wsolution);
    isoline.push_back(s.back());
  }

  std::string isofile("output/nurbs_solid_iso.txt");
  IO::writeXYdata(iso,isoline,isofile);
  std::system(std::string(python + "python/plot_xy.py " + isofile).c_str());
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

int main(int argc, char **argv)
{
  NurbsSolid solid;
  elbow(1,2,3,solid);
  // pipe(1,2,solid);
  simulation(solid);
  return 0;
}
