#include "splines/Algorithms.h"
#include "splines/BSpline.h"
#include "splines/Nurbs.h"
#include "splines/SplineModifiers.h"
#include "splines/utils/VectorOperations.h"
#include "splines/utils/Transformations.h"
#include "splines/utils/Converters.h"
#include "splines/utils/Quaternion.h"

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
using namespace vector_ops;

template<typename Curve>
void showCurve(Curve const &c)
{
  std::string const file = "output/curve.txt";
  spline_ops::writeToFile(c,file,10);
  std::system(std::string(python + "python/plot_curve.py " + file).c_str());
}

template<typename Surface>
void show(Surface const &s)
{
  std::string const file = "output/shell.txt";
  spline_ops::writeToFile(s,file,20,20);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());
}

void aplane(NurbsSurface &surface, double L = 1.0)
{
  surface.p = 1;
  surface.q = 1;
  surface.uknot = {0,0,1,1};
  surface.vknot = {0,0,1,1};
  surface.Q = {{0,0,0},{L,0,0}, {0,L,0}, {L,L,0}};
  surface.weights = {1,1,1,1};
}

void bplane(NurbsSurface &surface, double L = 1.0)
{
  surface.p = 1;
  surface.q = 1;
  surface.uknot = {0,0,1,1};
  surface.vknot = {0,0,1,1};
  surface.Q = {{0,0,0},{L,0,0}, {0,L,0}, {0,L,L}};
  surface.weights = {1,1,1,1};
}

void acylinder(NurbsSurface &surface, double R = 1.0)
{
  surface.p = 2;
  surface.q = 1;
  surface.uknot = {0,0,0,1,1,1};
  surface.vknot = {0,0,1,1};
  surface.Q = {{R,0,0},{R,R,0},{0,R,0},
               {R,0,R},{R,R,R},{0,R,R}};

  double const fact = std::sqrt(1.0/2.0);
  surface.weights = {1,fact,1,1,fact,1};
}

struct Plane
{
  std::vector<double> x;
  std::vector<double> n;
};

std::vector<size_t> pointsOnPlane(GeometricDofManager const &geom, Plane const &plane)
{
  auto on_plane = [&plane](auto const &p)
  {
    return algo::equal(dot(plane.n,p-plane.x),0.0);
  };

  std::vector<size_t> found;
  for (auto const id : geom.ids)
  {
    if (on_plane(geom.ctrlPoints[id])) found.push_back(id);
  }

  return found;
}

Eigen::SparseMatrix<double> assemblyMatrix(std::vector<size_t> const &sdof, size_t ndof, size_t ctrlPoints)
{
  auto const N = ndof*ctrlPoints;
  auto const M = ndof*sdof.size();
  Eigen::SparseMatrix<double> L(N,M);
  size_t k = 0;
  for (auto const &i :sdof)
  {
    for (size_t j = 0; j < ndof; j++)
    {
      L.insert(ndof*i+j,k++) = 1;
    }
  }

  return L;
}

template<typename Surface>
void run(Surface const &surf)
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
  // define weak forms
  double sigma = 400.0;
  double h = 1e-2;

  auto matl = [sigma](auto const &p)->Eigen::MatrixXd
  {
    Eigen::Matrix2d D(Eigen::Matrix2d::Identity());
    return sigma*D;
  };

  auto B = [](auto const &p)->SparseMatrix<double>
  {
    auto const grad   = p.mapper.grad(p.para);
    return grad.sparseView();
  };

  auto stiffness = [&B,&matl,h](auto const &p)->Eigen::MatrixXd
  {
    auto const d = matl(p);
    auto const b = B(p);

    auto const kel = h*b.transpose()*(d*b);
    return kel;
  };

  // set the up Simulation using the defined geometry
  GeometricDofManager geom;
  geom.addShape(&surf);

  // create mapper
  // // Compute stiffness
  size_t ndof = 1;
  size_t dof = ndof*geom.ctrlPoints.size();
  Eigen::MatrixXd K(dof,dof); K.setZero();
  Eigen::VectorXd f(dof); f.setZero();
  // integrate
  {
    ManifoldElementMapper mapper(surf);
    auto ip  = iga::integrationPoints(mapper,surf.p,surf.q,1);

    auto const L = assemblyMatrix(geom.idsForShape(&surf),ndof,geom.ctrlPoints.size());
    // integrate
    auto elu = iga::meshFromSpan(surf.uknot).size()-1;
    auto elv = iga::meshFromSpan(surf.vknot).size()-1;
    std::printf("Default mesh has {%lu, %lu} elements\n",elu,elv);

    size_t edof = ndof*surf.Q.size();
    Eigen::MatrixXd Kel(edof,edof);

    for (size_t j = 0; j < elv; j++)
    {
      for (size_t i = 0; i < elu; i++)
      {
        // std::printf("\nIntegrating element(%lu,%lu)\n",i,j);
        // std::cout << "updating element mesh..." << std::endl;
        mapper.updateElementMesh(i,j);
        // std::cout << "updating integration points..." << std::endl;
        for (auto &p : ip) p.update();
        // std::cout << "integrating weak forms..." << std::endl;
        Kel.setZero();
        quadrature::gauss(ip,Kel,stiffness);
        K += L*Kel*L.transpose();
      }
    }
  } 

  // apply boundary conditions
  std::vector<std::pair<size_t,double>> cdof;
  {
    Plane plane;
    plane.x = {0.0,0.0,0.0};
    plane.n = {0,-1,0};
    for (auto const &idx : pointsOnPlane(geom,plane))
    {
      cdof.push_back(std::make_pair(idx,300.0));   // fix all
    }

    plane.x = {0.0,1.0,0.0};
    plane.n = {0,1,0};
    for (auto const &idx : pointsOnPlane(geom,plane))
    {
      cdof.push_back(std::make_pair(idx,500.0));   // fix all
    }
  }

  size_t C = cdof.size();
  Eigen::MatrixXd Kc(C,dof); Kc.setZero();
  Eigen::VectorXd Rc(C); Rc.setZero();

  size_t k = 0;
  for (auto const &pr: cdof)
  {
    Kc(k,pr.first) = 1;
    Rc[k] = pr.second;
    k++;
  }

// #################################
// compose matrix 
// #################################
  Eigen::MatrixXd Kg(dof+C,dof+C); Kg.setZero();
  Eigen::VectorXd F(dof+C); F.setZero();
  F << f,Rc;
  Kg.block(0,0,dof,dof) = K;
  Kg.block(dof,0,C,dof) = Kc;
  Kg.block(0,dof,dof,C) = Kc.transpose();

  // IO::spy(K);
  // std::cout << "k = \n" << K << std::endl;
  // std::cout << "\nkg = \n" << Kg << std::endl;
  // std::cout << "\nkc = \n" << Kc << std::endl;
  // std::cout << "F = \n" << f.transpose() << std::endl;
  std::cout << "sum(F) = \n" << f.sum() << std::endl;
// #################################
// solve 
// #################################
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(Kg.sparseView());
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
// solution 
// #################################
  // std::cout << "**** Analytic solution = " << -0.58 << std::endl;

// #################################
// write solution to file and plot
// #################################
  // reshape solution
  // std::cout << "\n*** Solution vector: \n" << T << std::endl;
  NurbsSurface wsolution; IO::geometryWithSolution(surf,convert::to<std::vector<double>>(T),wsolution);
  show(wsolution);

  size_t N = 100;
  std::vector<double> v(N+1); std::iota(v.begin(),v.end(),0); v/= N;
  double u = 0.5;
  std::vector<double> isoline;
  for (auto x : v)
  {
    auto s = spline_ops::SurfacePoint(u,x,wsolution);
    isoline.push_back(s.back());
  }
  std::string isofile("output/nurbs_heat_iso.txt");
  IO::writeXYdata(v,isoline,isofile);
  std::system(std::string(python + "python/plot_xy.py " + isofile).c_str());
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

int main(int argc, char **argv)
{
  // NurbsSurface surf; aplane(surf,1.0);
  NurbsSurface surf; bplane(surf,1.0);
  show(surf);
  // NurbsCurve curve; bline(curve,1.0);
  spline_ops::elevate(2,0,surf);
  // spline_ops::elevate(1,1,surf);
  spline_ops::midpointRefinement(3,0,surf);
  spline_ops::midpointRefinement(1,1,surf);
  run(surf);
  return 0;
}
