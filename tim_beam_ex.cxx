#include "splines/Algorithms.h"
#include "splines/BSpline.h"
#include "splines/Nurbs.h"
#include "splines/SplineModifiers.h"
#include "splines/utils/VectorOperations.h"
#include "splines/utils/Transformations.h"
#include "splines/utils/Converters.h"
#include "splines/utils/Quaternion.h"

#include "iga/CurveElementMapper.h"
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

void arc(NurbsCurve &curve, double radius=1.0)
{
  using namespace vector_ops;

  int p = 2;
  std::vector<double> knot{0.0, 0.0, 0.0, 1, 1, 1};
  algo::normalizeKnot(knot);

  typename NurbsCurve::matrix cpts{
    {0.0, radius},
    { radius, radius },
    { radius , 0.0},
  };

  static double const fact{1.0/std::sqrt(2.0)};
  decltype(knot) weights
  {
    1.0, fact, 1.0
  };

  curve.p = p;
  curve.knot = knot;
  curve.weights = weights;
  curve.Q = cpts;

}

void aline(NurbsCurve &curve, double L = 1.0)
{
  curve.p = 1;
  curve.knot = {0,0,1,1};
  curve.Q = {{0,0},{L,0}};
  curve.weights = {1,1};
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

template<typename Curve>
void run(Curve const &curve, double factor)
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;

  // define weak forms
  double E = 1e6;
  double nu = 0.25;
  double h = 1.0/factor;
  double w = 1.0;
  double I = w*std::pow(h,3)/12.0;
  double A = h*w;
  double G = E/(2.0*(1.0+nu));

  auto Dmat = [E,nu,I,A,G](auto const &p)->Matrix2d
  {
    // double const kap = 10.0*(1.0+nu)/(12.0 + 11.0*nu);
    double const kap = 5.0/6.0;
    Matrix2d D(Matrix2d::Zero());
    D(0,0) = G*A*kap;
    D(1,1) = E*I;
    return D;
  };

  auto N = [](auto const &p)->MatrixXd
  {
    auto const shape = p.mapper.shape(p.para);
    auto const dof = shape.size();
    MatrixXd n(2,2*dof); n.setZero();
    for (int i = 0; i < dof; i ++)
    {
      n(0,i*2)   = shape[i];
      n(1,i*2+1) = shape[i];
    }

    return n;
  };

  auto B = [](auto const &p)->MatrixXd
  {
    auto const shape  = p.mapper.shape(p.para);
    auto const grad   = p.mapper.grad(p.para);
    auto const dof = shape.size();

    MatrixXd b(2,2*dof); b.setZero();
    for (int i = 0; i < dof; i++)
    {
      b(0,2*i) = grad(0,i);
      b(0,2*i+1) = shape(i);

      b(1,2*i+1) = grad(0,i);
    }
    return b;
  };

  auto stiffness = [&B,&Dmat](auto const &p)->Eigen::MatrixXd
  {
    auto const d = Dmat(p);
    auto const b = B(p);
    return b.transpose()*(d*b);
  };

  Eigen::Vector2d force; force << -1,0.0;
  auto load = [&N, &force](auto const &p)->Eigen::VectorXd
  {
    return N(p).transpose()*force;
    // return N(p).transpose()*(1.0 - p.para[0])*force;
  };

  // set the up Simulation using the defined geometry
  GeometricDofManager geom;
  geom.addShape(&curve);

  // create mapper
  // // Compute stiffness
  size_t ndof = 2;
  size_t dof = ndof*geom.ctrlPoints.size();
  Eigen::MatrixXd K(dof,dof); K.setZero();
  Eigen::VectorXd f(dof); f.setZero();
  // integrate
  {
    CurveElementMapper mapper(curve);
    auto ip = iga::integrationPoints(mapper,curve.p);

    auto const L = assemblyMatrix(geom.idsForShape(&curve),ndof,geom.ctrlPoints.size());

    // integrate
    auto elu = iga::meshFromSpan(curve.knot).size()-1;
    std::printf("Default mesh has {%lu} elements\n",elu);

    size_t edof = ndof*curve.Q.size();
    Eigen::MatrixXd Kel(edof,edof);
    Eigen::VectorXd fel(edof);

    for (size_t i = 0; i < elu; i++)
    {
      std::printf("\nIntegrating element(%lu)\n",i);
      // std::cout << "updating element mesh..." << std::endl;
      mapper.updateElementMesh(i);
      // std::cout << "updating integration points..." << std::endl;
      for (auto &p : ip) p.update();
      // std::cout << "integrating weak forms..." << std::endl;
      Kel.setZero(); fel.setZero();
      quadrature::gauss(ip,Kel,stiffness);
      quadrature::gauss(ip,fel,load);
      K += L*Kel*L.transpose();
      f += L*fel;
    }

  } 

  // apply boundary conditions
  std::vector<size_t> cdof;
  {
    Plane plane;
    plane.x = {0.5,0};
    plane.n = {1,0};
    for (auto const &idx : pointsOnPlane(geom,plane))
    {
      cdof.push_back(ndof*idx); // fix w
      // cdof.push_back(ndof*idx+1); // fix rot
    }

    plane.x = {0,0};
    plane.n = {-1,0};
    for (auto const &idx : pointsOnPlane(geom,plane))
    {
      cdof.push_back(ndof*idx+1); // fix rot
    }

    std::sort(cdof.begin(), cdof.end());
    auto last = std::unique(cdof.begin(),cdof.end());
    cdof.erase(last,cdof.end());
    // std::cout << "cdof = " << cdof << std::endl;
  }

  size_t C = cdof.size();
  Eigen::MatrixXd Kc(C,dof); Kc.setZero();
  Eigen::VectorXd Rc(C); Rc.setZero();

  size_t k = 0;
  for (auto const &i: cdof) Kc(k++,i) = 1;

// #################################
// compose matrix 
// #################################
  Eigen::MatrixXd Kg(dof+C,dof+C); Kg.setZero();
  Eigen::VectorXd F(dof+C); F.setZero();
  F << f,Rc;
  Kg.block(0,0,dof,dof) = K;
  Kg.block(dof,0,C,dof) = Kc;
  Kg.block(0,dof,dof,C) = Kc.transpose();

  // std::cout << "k = \n" << K << std::endl;
  // std::cout << "\nkg = \n" << Kg << std::endl;
  // std::cout << "\nkc = \n" << Kc << std::endl;
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

  Eigen::VectorXd const u = solver.solve(F).head(dof);
  if(solver.info()!=Eigen::ComputationInfo::Success)
  {
    std::cout << "Solve failed" << std::endl;
    return;
  }

// #################################
// solution 
// #################################
  std::cout << "**** Analytic solution = " << -0.58 << std::endl;

// #################################
// write solution to file and plot
// #################################
  // reshape solution
  auto const uv = Eigen::Map<const Eigen::MatrixXd>(u.data(),ndof,dof/ndof).transpose();
  std::cout << uv << std::endl;
  // compute disp mag and deform the mesh
  auto const sids = geom.idsForShape(&curve);
  Eigen::VectorXd vdisp(sids.size());

  double scale = E*std::pow(h,3.0);
  Curve deformed = curve; // make a copy
  for (size_t i = 0; i < sids.size(); i++)
  {
    auto const u = Eigen::RowVectorXd(uv.row(sids[i]));
    vdisp[i] = u[0];
    deformed.Q[i][1] += scale*u[0];
  }

  std::cout << "wdisp*" << scale << " = \n\n" << uv.col(0)*scale << std::endl;
  // plot deformed beam
  show(deformed);
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

template<typename Curve>
void show(Curve const &c)
{
  std::string const file = "output/curve.txt";
  spline_ops::writeToFile(c,file,20);
  std::system(std::string(python + "python/plot_curve.py " + file).c_str());
}

int main(int argc, char **argv)
{
  // NurbsCurve curve; arc(curve,1.0);
  NurbsCurve curve; aline(curve,0.5);
  spline_ops::elevate(3,curve);
  spline_ops::midpointRefinement(5,curve);
  run(curve,10.0);
  return 0;
}
