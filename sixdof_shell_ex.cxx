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

static Eigen::IOFormat clean(1, 0, ", ", "\n", "[", "]");

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
  spline_ops::writeToFile(s,file,10,10);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());
}

void incline(NurbsSurface &surface, double L = 1.0)
{
  surface.p = 1;
  surface.q = 1;
  surface.uknot = {0,0,1,1};
  surface.vknot = {0,0,1,1};
  surface.Q = {{0,0,0},{L,0,L},
              {0,L,0} ,{L,L,L}};
  surface.weights = {1,1,1,1};
}

void aplaneXY(NurbsSurface &surface, double L = 1.0)
{
  surface.p = 1;
  surface.q = 1;
  surface.uknot = {0,0,1,1};
  surface.vknot = {0,0,1,1};
  surface.Q = {{0,0,0},{L,0,0},
              {0,L,0} ,{L,L,0}};
  surface.weights = {1,1,1,1};
}

void aplaneXZ(NurbsSurface &surface, double L = 1.0)
{
  surface.p = 1;
  surface.q = 1;
  surface.uknot = {0,0,1,1};
  surface.vknot = {0,0,1,1};
  surface.Q = {{0,0,0},{L,0,0},
              {0,0,L} ,{L,0,L}};
  surface.weights = {1,1,1,1};
}

void aplaneYZ(NurbsSurface &surface, double L = 1.0)
{
  surface.p = 1;
  surface.q = 1;
  surface.uknot = {0,0,1,1};
  surface.vknot = {0,0,1,1};
  surface.Q = {{0,0,0},{0,L,0},
              {0,0,L} ,{0,L,L}};
  surface.weights = {1,1,1,1};
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
  size_t const ndof = 6;
  double const   h  = 0.01;
  double const   E = 70e9;
  double const  nu = 0.2;

  auto matl_membrane = [E,nu](auto const &p)->Eigen::Matrix3d
  {
    size_t const ik[] = {0,1,0};
    size_t const jl[] = {0,1,1};
    ManifoldElementMapper const * mapper =
      dynamic_cast<ManifoldElementMapper const *>(&p.mapper);

    auto const mu  = E/(2.0*(1+nu));
    // auto const lam = E*nu/((1+nu)*(1-2.0*nu));
    auto const lam = E*nu/((1+nu)*(1-nu));

    auto const A = mapper->contravariantMetricTensor(p.para);

    Eigen::Matrix3d D(Eigen::Matrix3d::Zero());
    for (size_t a = 0; a < 3; a++)
    {
      auto i = ik[a];
      auto j = jl[a];
      for (size_t b = 0; b < 3; b++)
      {
        auto k = ik[b];
        auto l = jl[b];
        D(a,b) = mu*(A(i,k)*A(j,l) + A(i,l)*A(j,k)) + lam*A(i,j)*A(k,l);
      }
    }
    return D;
  };

  auto matl_shear = [E,nu](auto const &p)->Eigen::Matrix2d
  {
    double const factor = 5.0/6.0;
    auto const mu  = E/(2.0*(1+nu));
    Eigen::Matrix2d D; D << mu,0,0,mu;
    return factor*D;
  };

  auto uid = [](VectorXd const X)
  {
    auto const dof = X.size();
    MatrixXd op(3,6*dof); op.setZero();
    for (int i = 0; i < dof; i ++)
    {
      for (size_t j = 0; j < 3; j++) op(j,6*i+j) = X[i];
    }
    return op;
  };

  auto rid = [](VectorXd const X)
  {
    auto const dof = X.size();
    MatrixXd op(3,6*dof); op.setZero();
    for (int i = 0; i < dof; i ++)
    {
      for (size_t j = 0; j < 3; j++)
        op(j,6*i+3+j) = X[i];
    }
    return op;
  };

  auto dot = [](auto const &a, auto const &b)
  {
    return a.transpose()*b;
  };

  // membrane strain operator
  auto Bm = [uid,dot](auto const &p)->SparseMatrix<double>
  {
    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&p.mapper);
    auto const grad = iga::CompactShapeFunctionDerivatives(p.para[0], p.para[1], mapper->manifold());
    auto const basis = mapper->covariantBasis(p.para);

    auto const dNd1 = uid(grad.row(0));
    auto const dNd2 = uid(grad.row(1));

    Vector3d const a1 = basis.col(0);
    Vector3d const a2 = basis.col(1);

    auto const dof = grad.cols();
    Eigen::MatrixXd b(3,6*dof); b.setZero();
    b.row(0) = dot(a1,dNd1);
    b.row(1) = dot(a2,dNd2);
    b.row(2) = dot(a1,dNd2) + dot(a2,dNd1);

    return b.sparseView();
  };

  // bending strain operator
  auto Bb = [rid,dot](auto const &p)->SparseMatrix<double>
  {
    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&p.mapper);
    auto const grad = iga::CompactShapeFunctionDerivatives(p.para[0], p.para[1], mapper->manifold());
    auto const basis = mapper->covariantBasis(p.para);

    auto const dNd1 = rid(grad.row(0));
    auto const dNd2 = rid(grad.row(1));

    Vector3d const a1 = basis.col(0);
    Vector3d const a2 = basis.col(1);

    auto const dof = grad.cols();
    Eigen::MatrixXd b(3,6*dof); b.setZero();
    b.row(0) = dot(a1,dNd1);
    b.row(1) = dot(a2,dNd2);
    b.row(2) = dot(a1,dNd2) + dot(a2,dNd1);
    return b.sparseView();
  };

  // transverse strain operator
  auto Bt = [uid,rid,dot](auto const &p)->SparseMatrix<double>
  {
    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&p.mapper);
    auto const shape = iga::CompactShapeFunctions(p.para[0], p.para[1], mapper->manifold());
    auto const grad  = iga::CompactShapeFunctionDerivatives(p.para[0], p.para[1], mapper->manifold());
    auto const basis = mapper->covariantBasis(p.para);

    auto const dNd1 = uid(grad.row(0));
    auto const dNd2 = uid(grad.row(1));
    auto const N    = rid(shape);

    Vector3d const a1 = basis.col(0);
    Vector3d const a2 = basis.col(1);
    Vector3d const n  = basis.col(2);

    auto const dof = grad.cols();
    Eigen::MatrixXd b(2,6*dof); b.setZero();

    b.row(0) = dot(n,dNd2) + dot(a2,N);
    b.row(1) = dot(n,dNd1) + dot(a1,N);

    return b.sparseView();
  };

  auto membrane_stiffness = [&Bm,&matl_membrane](auto const &p)->Eigen::MatrixXd
  {
    auto const d = matl_membrane(p);
    auto const b = Bm(p);

    return b.transpose()*d*b;
  };

  auto bending_stiffness = [&Bb,&matl_membrane](auto const &p)->Eigen::MatrixXd
  {
    auto const d = matl_membrane(p);
    auto const b = Bb(p);

    return b.transpose()*d*b;
  };

  auto transverse_stiffness = [&Bt,&matl_shear](auto const &p)->Eigen::MatrixXd
  {
    auto const d = matl_shear(p);
    auto const b = Bt(p);
    return b.transpose()*d*b;
  };

  auto constraint = [dot,rid](const auto &p)->MatrixXd
  {
    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&p.mapper);
    auto const shape = p.mapper.shape(p.para);
    auto const n = mapper->covariantBasis(p.para).col(2);
    auto const x = dot(n,rid(shape));
    return 1e8*x.transpose()*x;
  };

  auto stiffness = [
    &membrane_stiffness,
    &bending_stiffness,
    &transverse_stiffness,
    &constraint,h](auto const &p)->Eigen::MatrixXd
  {
    auto kc = constraint(p); 
    auto km = membrane_stiffness(p);
    // std::cout << "km = \n" << km.format(clean) << std::endl;
    auto kb = bending_stiffness(p);
    // std::cout << "\nkb = \n" << kb.format(clean) << std::endl;
    auto kt = transverse_stiffness(p);
    // std::cout << "\nkt = \n" << kt.format(clean) << std::endl;
    return h*(km+kt) + std::pow(h,3)/12.0*kb + kc;
  };

  double const pressure = 25e3;
  auto load = [&uid, pressure](auto const &p)->Eigen::VectorXd
  {
    auto const mapper = dynamic_cast<ManifoldElementMapper const *>(&p.mapper);
    auto const n = mapper->covariantBasis(p.para).col(2);
    auto const shape = p.mapper.shape(p.para);
    return -pressure*uid(shape).transpose()*n;
  };

  // set the up Simulation using the defined geometry
  GeometricDofManager geom;
  geom.addShape(&surf);

  // create mapper
  // // Compute stiffness
  size_t dof = ndof*geom.ctrlPoints.size();
  Eigen::MatrixXd K(dof,dof); K.setZero();
  Eigen::VectorXd f(dof); f.setZero();
  // integrate
  {
    ManifoldElementMapper mapper(surf);
    auto ip  = iga::integrationPoints(mapper,surf.p,surf.q);

    auto const L = assemblyMatrix(geom.idsForShape(&surf),ndof,geom.ctrlPoints.size());
    // integrate
    auto elu = iga::meshFromSpan(surf.uknot).size()-1;
    auto elv = iga::meshFromSpan(surf.vknot).size()-1;
    std::printf("Default mesh has {%lu, %lu} elements\n",elu,elv);

    size_t edof = ndof*surf.Q.size();
    Eigen::MatrixXd Kel(edof,edof);
    Eigen::VectorXd fel(edof);

    for (size_t j = 0; j < elv; j++)
    {
      for (size_t i = 0; i < elu; i++)
      {
        std::printf("\nIntegrating element(%lu,%lu)\n",i,j);
        // std::cout << "updating element mesh..." << std::endl;
        mapper.updateElementMesh(i,j);
        // std::cout << "updating integration points..." << std::endl;
        for (auto &p : ip) p.update();
        // std::cout << "integrating weak forms..." << std::endl;
        Kel.setZero();
        quadrature::gauss(ip,Kel,stiffness);
        K += L*Kel*L.transpose();
        fel.setZero();
        quadrature::gauss(ip,fel,load);
        f += L*fel;
      }
    }
  } 

  // apply boundary conditions
  std::vector<size_t> cdof;
  {
    Plane plane;
    plane.x = {0.0,0.0,0.0};
    plane.n = {-1,0,0};
    for (auto const &idx : pointsOnPlane(geom,plane))
    {
      for (size_t i = 0; i < ndof; i++)
        cdof.push_back(ndof*idx+i);   // fix all
    }

    std::sort(cdof.begin(), cdof.end());
    auto last = std::unique(cdof.begin(),cdof.end());
    cdof.erase(last,cdof.end());
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

  std::cout << "stiffness: \n" << std::endl;
  std::cout << K.format(clean) << std::endl;
  // std::cout << "load: \n" << std::endl;
  // std::cout << f.format(clean) << std::endl;
  // IO::spy(Kg);
  // std::cout << "k = \n" << K << std::endl;
  // std::cout << "\nkg = \n" << Kg << std::endl;
  // std::cout << "\nkc = \n" << Kc << std::endl;
  // std::cout << "F = \n" << f.transpose() << std::endl;
  // std::cout << "sum(F) = \n" << f.sum() << std::endl;
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
  // std::cout << "**** Analytic solution = " << -0.58 << std::endl;

// #################################
// write solution to file and plot
// #################################
  // reshape solution
  MatrixXd const uvw = Eigen::Map<const Eigen::MatrixXd>(u.data(),ndof,dof/ndof).transpose();
  std::cout << "\n*** Solution vector: \n" << uvw << std::endl;
  // compute disp mag and deform the mesh
  auto const sids = geom.idsForShape(&surf);
  Eigen::VectorXd  umag(sids.size());
  Eigen::VectorXd  w(sids.size());

  // double mx = std::max(uvw.maxCoeff(), std::abs(uvw.minCoeff()));
  // double scale = mx > 0 ? 1.0/mx : 1.0;
  double scale = 1.0;
  std::cout << "Scale set to = " << scale << std::endl; 

  Surface deformed = surf; // make a copy
  for (size_t i = 0; i < sids.size(); i++)
  {
    auto const u = Eigen::RowVectorXd(uvw.row(sids[i]));
    umag[i]  = u.head(3).norm();
    w[i]     = u[2];
    deformed.Q[i][0] += scale*u[0]; // add displacement to coordinates
    deformed.Q[i][1] += scale*u[1];
    deformed.Q[i][2] += scale*u[2];
  }

  auto range = [](auto const &x)->std::pair<double,double>
  {
    return std::make_pair(x.minCoeff(), x.maxCoeff());
  };

  std::cout << "umag = " << range(umag) << std::endl;
  std::cout << "w    = " << range(w) << std::endl;
  // NurbsCurve iso;
  // spline_ops::getIsoCurve(0,1,deformed,iso);
  // std::cout << iso << std::endl;
  // std::cout << "wdisp*" << scale << " = \n\n" << uv.col(1)*scale << std::endl;
  // plot deformed beam
  // showCurve(iso);
  show(deformed);
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

int main(int argc, char **argv)
{
  NurbsSurface surf;
  // incline(surf,1.0/sqrt(2.0));
  aplaneXY(surf,1.0);
  // show(surf);
  // NurbsCurve curve; bline(curve,1.0);
  // spline_ops::elevate(2,0,surf);
  // spline_ops::elevate(2,1,surf);
  // spline_ops::midpointRefinement(2,0,surf);
  // spline_ops::midpointRefinement(2,1,surf);
  run(surf);
  return 0;
}
