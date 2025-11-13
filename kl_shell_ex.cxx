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



void aplane(NurbsSurface &surface, double L = 1.0)
{
  surface.p = 2;
  surface.q = 2;
  surface.uknot = {0,0,0,1,1,1};
  surface.vknot = {0,0,0,1,1,1};
  surface.Q = {{0,0,0} ,{L/2,0,0}  ,{L,0,0},
              {0,L/2,0},{L/2,L/2,0},{L,L/2,0},
              {0,L,0}  ,{L/2,L,0}  ,{L,L,0}};
  surface.weights = {1,1,1,1,1,1,1,1,1};
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
  double nu = 0.2;
  double h  = 0.01;

  auto matl = [E,nu](auto const &p)->Eigen::Matrix3d
  {
    size_t const ik[] = {0,1,0};
    size_t const jl[] = {0,1,1};

    auto const mu  = E/(2.0*(1+nu));
    auto const lam = E*nu/((1+nu)*(1-2.0*nu));

    auto const A = p.mapper.contravariantMetricTensor(p.para);

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

  auto id = [](VectorXd const X)
  {
    auto const dof = X.size();
    MatrixXd op(3,3*dof); op.setZero();
    for (int i = 0; i < dof; i ++)
    {
      for (size_t j = 0; j < 3; j++) op(j,3*i+j) = X[i];
    }
    return op;
  };

  auto dot = [](auto const &a, auto const &b)
  {
    return a.transpose()*b;
  };

  auto N = [id](auto const &p)->MatrixXd
  {
    auto const shape = p.mapper.shape(p.para);
    return id(shape);
  };

  // membrane strain operator
  auto Bm = [id, dot](auto const &p)->SparseMatrix<double>
  {
    auto const grad = iga::ShapeFunctionDerivatives(p.para[0], p.para[1], p.mapper.manifold());
    auto const basis = p.mapper.covariantBasis(p.para);

    auto const dNd1 = id(grad.row(0));
    auto const dNd2 = id(grad.row(1));

    Vector3d const a1 = basis.col(0);
    Vector3d const a2 = basis.col(1);

    auto const dof = grad.cols();
    Eigen::MatrixXd b(3,3*dof);

    b.row(0) = dot(a1,dNd1);
    b.row(1) = dot(a2,dNd2);
    b.row(2) = dot(a2,dNd1) + dot(a1,dNd2);
    return b.sparseView();
  };

  // bending strain operator
  auto Bb = [id, dot](auto const &p)->SparseMatrix<double>
  {
    auto const grad = iga::ShapeFunctionDerivatives(p.para[0], p.para[1], p.mapper.manifold());
    auto const basis = p.mapper.covariantBasis(p.para);

    MatrixXd const dNdu = id(grad.row(0));
    MatrixXd const dNdv = id(grad.row(1));

    MatrixXd    const Q     = convert::to<MatrixXd>(p.mapper.manifold().Q);
    RowVectorXd const dNduu = iga::ShapeFunctionDerivative2(p.para[0],p.para[1],0,p.mapper.manifold());
    RowVectorXd const dNdvv = iga::ShapeFunctionDerivative2(p.para[0],p.para[1],1,p.mapper.manifold());
    RowVectorXd const dNduv = iga::ShapeFunctionMixedDerivative(p.para[0],p.para[1],p.mapper.manifold());

    Vector3d const a1 = basis.col(0);
    Vector3d const a2 = basis.col(1);
    Vector3d const n  = basis.col(2);

    double const mag = std::sqrt(a1.dot(a1)*a2.dot(a2) - std::pow(a1.dot(a2), 2));

    Vector3d const da1d1 = dNduu*Q; 
    Vector3d const da2d2 = dNdvv*Q; 
    Vector3d const da1d2 = dNduv*Q; 

    auto const dof = grad.cols();
    Eigen::MatrixXd b(3,3*dof);

    b.row(0) = -dot(n,id(dNduu))
             + (dot(da1d1.cross(a2),dNdu) + dot(a1.cross(da1d1),dNdv))/mag
             + (dot(a2.cross(n),dNdu) + dot(n.cross(a1),dNdv))*(n.dot(da1d1)/mag);
    b.row(1) = -dot(n,id(dNdvv))
             + (dot(da2d2.cross(a2),dNdu) + dot(a1.cross(da2d2),dNdv))/mag
             + (dot(a2.cross(n),dNdu) + dot(n.cross(a1),dNdv))*(n.dot(da2d2)/mag);
    b.row(2) = -2.0*dot(n,id(dNduv))
             + 2.0*(dot(da1d2.cross(a2),dNdu) + dot(a1.cross(da1d2),dNdv))/mag
             + 2.0*(dot(a2.cross(n),dNdu) + dot(n.cross(a1),dNdv))*(n.dot(da1d2)/mag);
    return b.sparseView();
  };

  // auto transformation = [](auto const &p, size_t ndof)->Eigen::SparseMatrix<double>
  // {
  //   MatrixXd R(MatrixXd::Zero(5,5));
  //   R.block(0,0,3,3) = p.mapper.localTransformation(p.para);
  //   R(3,3) = 1.0;
  //   R(4,4) = 1.0;

  //   MatrixXd th(MatrixXd::Zero(ndof,ndof));

  //   for (size_t i = 0; i < ndof/5; i++)
  //   {
  //     th.block(5*i,5*i,5,5) = R;
  //   }

  //   return th.sparseView();
  // };
  // 


  auto membrane_stiffness = [&Bm,&matl,h](auto const &p)->Eigen::MatrixXd
  {
    auto const d = matl(p);
    auto const b = Bm(p);

    return h*(b.transpose()*d*b);
  };

  auto bending_stiffness = [&Bb,&matl,h](auto const &p)->Eigen::MatrixXd
  {
    auto const d = matl(p);
    auto const b = Bb(p);

    return (std::pow(h,3)/12.0)*(b.transpose()*d*b);
  };

  double const pressure = 25e3;
  auto load = [&N, pressure](auto const &p)->Eigen::VectorXd
  {
    auto const n = p.mapper.covariantBasis(p.para).col(2);
    return -pressure*N(p).transpose()*n;
  };

  // set the up Simulation using the defined geometry
  GeometricDofManager geom;
  geom.addShape(&surf);

  // create mapper
  // // Compute stiffness
  size_t ndof = 3;
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
        quadrature::gauss(ip,Kel,membrane_stiffness);
        K += L*Kel*L.transpose();

        Kel.setZero();
        quadrature::gauss(ip,Kel,bending_stiffness);
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

    // something way fucky == no rotation at the edge
    for (size_t j = 0; j < surf.vknot.size()-surf.q-1; j++)
    {
      auto idx = surf.qid(1,j);
      cdof.push_back(ndof*idx+2);   // fix all
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

  // std::cout << "stiffness: \n" << std::endl;
  // std::cout << K.format(clean) << std::endl;
  // std::cout << "load: \n" << std::endl;
  // std::cout << f.format(clean) << std::endl;
  // IO::spy(Kg);
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
  NurbsSurface surf; aplane(surf,1.0);
  // show(surf);
  // NurbsCurve curve; bline(curve,1.0);
  spline_ops::elevate(2,0,surf);
  spline_ops::elevate(2,1,surf);
  spline_ops::midpointRefinement(4,0,surf);
  spline_ops::midpointRefinement(1,1,surf);
  run(surf);
  return 0;
}
