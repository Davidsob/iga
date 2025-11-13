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
  spline_ops::writeToFile(s,file,10,10);
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
  double E = 70e9;
  double nu = 0.33;
  double h = 1e-2;
  double kap = 5.0/6.0;

  auto matl = [E,nu,kap](auto const &p)->Eigen::MatrixXd
  {
    Eigen::MatrixXd D(6,6); D.setZero();
    double const a = E/(1.0-nu*nu);
    double const b = 0.5*(1.0-nu);

    D(0,0) = 1.0;
    D(0,1) = nu;

    D(1,0) = nu;
    D(1,1) = 1.0;


    D(3,3) = b;
    D(4,4) = kap*b;
    D(5,5) = kap*b;

    return a*D;
  };

  auto N = [h](auto const &p)->MatrixXd
  {
    auto const shape = p.mapper.shape(p.para);
    auto const dof = shape.size();
    MatrixXd n(3,5*dof); n.setZero();
    for (int i = 0; i < dof; i ++)
    {
      n(0,5*i)   = shape[i];
      n(1,5*i+1) = shape[i];
      n(2,5*i+2) = shape[i];
    }

    return n;
  };

  auto dudx = [h](auto const &p)->MatrixXd
  {
    auto const shape  = p.mapper.shape(p.para);
    auto const grad   = p.mapper.grad(p.para);

    auto const dof  = shape.size();
    auto const ndof = 5;

    auto const tau = p.mapper.tangents(p.para);
    auto const &g1 = -0.5*h*tau.row(1);
    auto const &g2 =  0.5*h*tau.row(0);
    auto const x   = p.gauss[2];

    MatrixXd du(9,ndof*dof); du.setZero();

    for (int i = 0; i < dof; i++)
    {
      // du
      du(0,ndof*i)   = grad(0,i);
      du(0,ndof*i+3) = grad(0,i)*x*g1[0];
      du(0,ndof*i+4) = grad(0,i)*x*g2[0];

      du(1,ndof*i)   = grad(1,i);
      du(1,ndof*i+3) = grad(1,i)*x*g1[0];
      du(1,ndof*i+4) = grad(1,i)*x*g2[0];

      du(2,ndof*i+3) = shape[i]*g1[0];
      du(2,ndof*i+4) = shape[i]*g2[0];

      // dv
      du(3,ndof*i+1) = grad(0,i);
      du(3,ndof*i+3) = grad(0,i)*x*g1[1];
      du(3,ndof*i+4) = grad(0,i)*x*g2[1];

      du(4,ndof*i+1) = grad(1,i);
      du(4,ndof*i+3) = grad(1,i)*x*g1[1];
      du(4,ndof*i+4) = grad(1,i)*x*g2[1];

      du(5,ndof*i+3) = shape[i]*g1[1];
      du(5,ndof*i+4) = shape[i]*g2[1];

      // dw
      du(6,ndof*i+2) = grad(0,i);
      du(6,ndof*i+3) = grad(0,i)*x*g1[2];
      du(6,ndof*i+4) = grad(0,i)*x*g2[2];

      du(7,ndof*i+2) = grad(1,i);
      du(7,ndof*i+3) = grad(1,i)*x*g1[2];
      du(7,ndof*i+4) = grad(1,i)*x*g2[2];

      du(8,ndof*i+3) = shape[i]*g1[2];
      du(8,ndof*i+4) = shape[i]*g2[2];
    }
    return du;
  };

  auto B = [&dudx](auto const &p)->SparseMatrix<double>
  {
    auto const du = dudx(p);

    MatrixXd b(6,du.cols()); b.setZero();
    b.row(0) = du.row(0);
    b.row(1) = du.row(4);
    b.row(2) = du.row(8);
    b.row(3) = du.row(1) + du.row(3); // dudy + dvdx
    b.row(4) = du.row(5) + du.row(7); // dvdz + dwdy
    b.row(5) = du.row(2) + du.row(6); // dudz + dwdx

    return b.sparseView();
  };

  auto transformation = [](auto const &p, size_t ndof)->Eigen::SparseMatrix<double>
  {
    MatrixXd R(MatrixXd::Zero(5,5));
    R.block(0,0,3,3) = p.mapper.localTransformation(p.para);
    R(3,3) = 1.0;
    R(4,4) = 1.0;

    MatrixXd th(MatrixXd::Zero(ndof,ndof));

    for (size_t i = 0; i < ndof/5; i++)
    {
      th.block(5*i,5*i,5,5) = R;
    }

    return th.sparseView();
  };

  auto stiffness = [&B,&matl,&transformation,h](auto const &p)->Eigen::MatrixXd
  {
    auto const d = matl(p);
    auto const b = B(p);

    auto const kel = h*b.transpose()*(d*b);
    return kel;

    // auto const R = transformation(p,kel.rows());
    // return R.transpose()*kel*R;
  };

  double pressure = 25e7;
  auto load = [&N, pressure](auto const &p)->Eigen::VectorXd
  {
    auto const n = p.mapper.normal(p.para);
    return -pressure*N(p).transpose()*n;
  };

  // set the up Simulation using the defined geometry
  GeometricDofManager geom;
  geom.addShape(&surf);

  // create mapper
  // // Compute stiffness
  size_t ndof = 5;
  size_t dof = ndof*geom.ctrlPoints.size();
  Eigen::MatrixXd K(dof,dof); K.setZero();
  Eigen::VectorXd f(dof); f.setZero();
  // integrate
  {
    ManifoldElementMapper mapper(surf);
    auto ip  = iga::integrationPoints(mapper,surf.p,surf.q,1);
    auto ips = iga::integrationPoints(mapper,surf.p,surf.q);

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
        for (auto &p : ips) p.update();
        quadrature::gauss(ips,fel,load);
        f += L*fel;
      }
    }
  } 


  // apply boundary conditions
  std::vector<size_t> cdof;
  {
    Plane plane;
    plane.x = {0.0,0.0,0.0};
    plane.n = {1,0,0};
    for (auto const &idx : pointsOnPlane(geom,plane))
    {
      for (size_t i = 0; i < ndof; i++)
        cdof.push_back(ndof*idx+i);   // fix all
    }

    std::sort(cdof.begin(), cdof.end());
    auto last = std::unique(cdof.begin(),cdof.end());
    cdof.erase(last,cdof.end());
    std::cout << cdof << std::endl;
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
  auto const uvw = Eigen::Map<const Eigen::MatrixXd>(u.data(),ndof,dof/ndof).transpose();
  std::cout << "\n*** Solution vector: \n" << uvw << std::endl;
  // compute disp mag and deform the mesh
  auto const sids = geom.idsForShape(&surf);
  Eigen::VectorXd  umag(sids.size());

  double scale = 1.0;
  Surface deformed = surf; // make a copy
  for (size_t i = 0; i < sids.size(); i++)
  {
    auto const u = Eigen::RowVectorXd(uvw.row(sids[i]));
    umag[i]  = u.head(3).norm();
    deformed.Q[i][0] += scale*u[0]; // add displacement to coordinates
    deformed.Q[i][1] += scale*u[1];
    deformed.Q[i][2] += scale*u[2];
  }

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
  NurbsSurface surf; aplane(surf,1.0);
  // show(surf);
  // NurbsCurve curve; bline(curve,1.0);
  spline_ops::elevate(2,0,surf);
  // spline_ops::elevate(2,1,surf);
  // spline_ops::midpointRefinement(2,0,surf);
  // spline_ops::midpointRefinement(2,1,surf);
  run(surf);
  return 0;
}
