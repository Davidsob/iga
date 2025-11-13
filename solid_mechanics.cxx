#include "splines/Algorithms.h"
#include "splines/BSpline.h"
#include "splines/Nurbs.h"
#include "splines/SplineModifiers.h"
#include "splines/utils/VectorOperations.h"
#include "splines/utils/Transformations.h"
#include "splines/utils/Converters.h"
#include "splines/utils/Quaternion.h"

#include "iga/GeometricDofManager.h"
#include "iga/IntegrationPoints.h"
#include "iga/IgaIO.h"
#include "iga/ManifoldElementMapper.h"
#include "iga/ParametricMesh.h"
#include "iga/ShapeFunctions.h"
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
#include <Eigen/SparseCholesky>

#include <random>

using namespace Eigen;

static const Eigen::IOFormat CleanFormat(2, 0, ", ", "\n", "[", "]");

void beam(NurbsSolid &solid, double width=1.0, double height=1.0, double depth=1.0)
{
  using namespace vector_ops;

  int p = 1;
  int q = 1;
  int r = 1;
  std::vector<double> knot{0.0, 0.0, 1.0, 1.0};
  algo::normalizeKnot(knot);

  typename NurbsCurve::matrix cpts{
    {0        , 0.0   , 0.0},
    {width    , 0.0   , 0.0},
    {0        , height, 0.0},
    {width    , height, 0.0},

    {0        , 0.0   , depth},
    {width    , 0.0   , depth},
    {0        , height, depth},
    {width    , height, depth},

  };

  solid.p = p;
  solid.q = q;
  solid.r = r;
  solid.uknot = knot;
  solid.vknot = knot;
  solid.wknot = knot;
  solid.weights = std::vector<double>(cpts.size(),1.0);
  solid.Q = cpts;
}

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

struct Plane
{
  std::vector<double> x;
  std::vector<double> n;
};

template<typename Shape>
std::vector<size_t> pointsOnPlaneForShape(Shape const &shape, Plane const &plane)
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
    L.insert(ndof*i,k++) = 1;
    L.insert(ndof*i+1,k++) = 1;
    L.insert(ndof*i+2,k++) = 1;
  }

  return L;
}

template<typename Solid>
void simulation(Solid const &solid)
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
// #################################
// Operators 
// #################################
  auto Bop = [](auto const &p)->Eigen::SparseMatrix<double>
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
      b(5,i*3+2) = grad(0,i);
    }

    return b.sparseView();
  };

  auto Nop = [](auto const &p)->Eigen::SparseMatrix<double>
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
    return N.sparseView();
  };

// #################################
// Material 
// #################################
  double const E = 70e9;
  double const nu = 0.29;
  // double const rho = 2700.0;
  auto matl = [E,nu](auto const &p)->Eigen::MatrixXd
  {
    Eigen::MatrixXd D(6,6); D.setZero();
    double const a = E/((1.0+nu)*(1.0-2.0*nu));
    double const b = 1.0-nu;
    double const c = 0.5*(1.0-2.0*nu);

    D(0,0) = b;
    D(0,1) = nu;
    D(0,2) = nu;

    D(1,0) = nu;
    D(1,1) = b;
    D(1,2) = nu;

    D(2,0) = nu;
    D(2,1) = nu;
    D(2,2) = b;

    D(3,3) = c;
    D(4,4) = c;
    D(5,5) = c;

    return a*D;
  };

// #################################
// Weak forms 
// #################################
  // Eigen::Vector3d accel; accel << 0.0, 0.0, -9.81;
  // auto gravity = [&Nop, &rho, &accel](auto const &p)->Eigen::VectorXd
  // {
  //   auto const n = Nop(p);
  //   return rho*(n.transpose()*accel);
  // };

  Eigen::Vector3d tapp; tapp << 0.0, 0.0, -1e8;
  auto traction = [&Nop, &tapp](auto const &p)->Eigen::VectorXd
  {
    auto const n = Nop(p);
    return n.transpose()*tapp;
  };

  auto stiffness = [&Bop,&matl](auto const &p)->Eigen::MatrixXd
  {
    Eigen::MatrixXd const b = Bop(p);
    Eigen::MatrixXd const d = matl(p);
    return b.transpose()*(d*b);
  };

  // set the up Simulation using the defined geometry
  GeometricDofManager geom;
  geom.addShape(&solid);

  // Compute stiffness
  size_t ndof = 3;
  size_t dof = ndof*geom.ctrlPoints.size();
  Eigen::MatrixXd K(dof,dof); K.setZero();
  Eigen::VectorXd f(dof); f.setZero();

  {
    size_t edof = ndof*solid.Q.size();
    Eigen::MatrixXd Kel(edof,edof);
    Eigen::VectorXd fel(edof);
    // create mapper for computation of volumetric terms
    SolidElementMapper mapper(solid);
    auto ip = iga::integrationPoints(mapper,solid.p,solid.q,solid.r);
    // get assembly matrix for shape
    auto const L = assemblyMatrix(geom.idsForShape(&solid),ndof,geom.ctrlPoints.size());
    // integrate  the volume terms
    auto elu = iga::meshFromSpan(solid.uknot).size()-1;
    auto elv = iga::meshFromSpan(solid.vknot).size()-1;
    auto elw = iga::meshFromSpan(solid.wknot).size()-1;
    std::printf("Default mesh has {%lu, %lu, %lu} elements\n",elu,elv,elw);
    for (size_t k = 0; k < elw; k++)
    {
      for (size_t j = 0; j < elv; j++)
      {
        for (size_t i = 0; i < elu; i++)
        {
          std::printf("\nIntegrating element(%lu,%lu,%lu)\n",i,j,k);
          std::cout << "updating element mesh..." << std::endl;
          mapper.updateElementMesh(i,j,k);
          std::cout << "updating integration points..." << std::endl;
          for (auto &p : ip) p.update();
          std::cout << "integrating weak forms..." << std::endl;
          Kel.setZero(); fel.setZero();
          quadrature::gauss(ip,Kel,stiffness);
          // quadrature::gauss(ip,fel,gravity);
          K += L*Kel*L.transpose();
          // f += L*fel;
        }
      }
    }
  }

  // integrate boundary terms
  {
    NurbsSurface manifold; spline_ops::getIsoSurface(1.0,2,solid,manifold);

    size_t edof = ndof*manifold.Q.size();
    Eigen::VectorXd fel(edof);

    ManifoldElementMapper mapper(manifold);
    geom.addShape(&manifold);
    // get some integration points for the manifold
    auto ip = iga::integrationPoints(mapper,manifold.p,manifold.q);
    // get assembly matrix for shape
    auto const L = assemblyMatrix(geom.idsForShape(&manifold),ndof,geom.ctrlPoints.size());
    // compute the surface integral
    auto elu = iga::meshFromSpan(manifold.uknot).size()-1;
    auto elv = iga::meshFromSpan(manifold.vknot).size()-1;
    std::printf("Default boundary mesh has {%lu, %lu} elements\n",elu,elv);
    
    for (size_t j = 0; j < elv; j++)
    {
      for (size_t i = 0; i < elu; i++)
      {
        std::printf("\nIntegrating boundary element(%lu,%lu)\n",i,j);
        std::cout << "updating boundary element mesh..." << std::endl;
        mapper.updateElementMesh(i,j);
        std::cout << "updating boundary integration points..." << std::endl;
        for (auto &p : ip) p.update();
        std::cout << "integrating boundary weak forms..." << std::endl;
        fel.setZero();
        quadrature::gauss(ip,fel,traction);
        f += L*fel;
      }
    }
    std::cout << "total load applied: " << f.sum() << std::endl;
  }

// #################################
// Boundary conditions
// #################################
  // compute constraints
  std::vector<size_t> cdof;
  // fix at the top
  Plane plane;
  plane.x = {0,0,4};
  plane.n = {0,0,1};
  // x- plane
  auto const fixed = pointsOnPlane(geom,plane);
  for (auto const &idx : fixed)
  {
    cdof.push_back(ndof*idx);
    cdof.push_back(ndof*idx+1);
    cdof.push_back(ndof*idx+2);
  }
  std::printf("Applying boundary conditions for fixed surface to %lu cpts\n",fixed.size());

  size_t C = cdof.size();
  Eigen::MatrixXd Kc(C,dof); Kc.setZero();
  Eigen::VectorXd Rc(C); Rc.setZero();

  size_t k = 0;
  for (auto const &i: cdof) Kc(k++,i) = 1;
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
  auto const uvw = Eigen::Map<const Eigen::MatrixXd>(u.data(),ndof,dof/ndof).transpose();
  // compute disp mag and deform the mesh
  auto const sids = geom.idsForShape(&solid);
  Eigen::VectorXd  umag(sids.size());
  Eigen::VectorXd wdisp(sids.size());

  double scale = 10;
  Solid deformed = solid;
  for (size_t i = 0; i < sids.size(); i++)
  {
    auto const x = convert::to<std::vector<double>>(Eigen::RowVectorXd(uvw.row(sids[i])));
    umag[i]  = norm(x);
    wdisp[i] = x[2];
    deformed.Q[i] += scale*x;
  }

  std::string file("output/nurbs_stress.txt");
  IO::writeSolutionToFile(deformed,wdisp,file,16,4,32);
  std::system(std::string(python + "python/plot_surface.py " + file).c_str());

  // check along v=w=const
  NurbsSolid wsolution; IO::geometryWithSolution(solid,wdisp,wsolution);
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

template<typename Solid>
void SolidTest(Solid const &solid)
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  using namespace vector_ops;

  int N = 5;
  double u,v,w,du,dv,dw;
  u = 1.0; v = 0.0; w = 0.0;
  du = 1.0/N, dv = 1.0/N, dw = 1.0/N;
  std::vector<std::vector<double>> p,pu,pv,pw;
  Eigen::MatrixXd const Q = convert::to<Eigen::MatrixXd>(solid.Q);
  for (int k = 0; k <= N; k++)
  {
    v = 0.0;
    for (int j = 0; j <= N; j++)
    {
      u = 0.0;
      for (int i = 0; i <= N; i++)
      {
        p.push_back(spline_ops::SolidPoint(u,v,w,solid));
        Eigen::MatrixXd const dN = iga::ShapeFunctionDerivatives(u,v,w,solid);
        Eigen::MatrixXd const grad = dN*Q;
        pu.push_back(convert::to<std::vector<double>>(Eigen::RowVectorXd(grad.row(0))));
        pv.push_back(convert::to<std::vector<double>>(Eigen::RowVectorXd(grad.row(1))));
        pw.push_back(convert::to<std::vector<double>>(Eigen::RowVectorXd(grad.row(2))));
        // pu.push_back(spline_ops::SolidDerivative(u,v,w,1,0,solid));
        // pv.push_back(spline_ops::SolidDerivative(u,v,w,1,1,solid));
        // pw.push_back(spline_ops::SolidDerivative(u,v,w,1,2,solid));
        u += du;
      }
      v += dv;
    }
    w += dw;
  }

  std::string file("output/nurbs_solid.txt");
  std::string uvec_file("output/nurbs_solid_du.txt");
  std::string vvec_file("output/nurbs_solid_dv.txt");
  std::string wvec_file("output/nurbs_solid_dw.txt");
  spline_ops::writeVectorData(p,pu, uvec_file, true, 0.2);
  spline_ops::writeVectorData(p,pv, vvec_file, true, 0.2);
  spline_ops::writeVectorData(p,pw, wvec_file, true, 0.2);
  spline_ops::writeToFile(solid,file,1,1,1);

  std::string cmd = python + "python/plot_surface.py ";
  cmd += file + " ";
  cmd += uvec_file + " ";
  cmd += vvec_file + " ";
  cmd += wvec_file;

  std::system(cmd.c_str());
}

template<typename Solid>
void testSurfaceIntegration(Solid const &solid)
{
  // get end of the pipe and integrate the surface area
  NurbsSurface manifold; spline_ops::getIsoSurface(0.0,1,solid,manifold);
  std::string const file("output/nurbs_surface.txt");
  spline_ops::writeToFile(manifold,file,20,20);

  // strategy put into an element mapper
  // get integration points appropriate for manifold surface
  // integrate
 
  // use a manifold element mapper 
  ManifoldElementMapper mapper(manifold);
  std::string cmd = python + "python/plot_surface.py ";
  cmd += file;
  std::system(cmd.c_str());
  // get some integration points for the manifold
  auto ip = iga::integrationPoints(mapper,manifold.p,manifold.q);
  // compute the surface integral
  auto elu = iga::meshFromSpan(manifold.uknot).size()-1;
  auto elv = iga::meshFromSpan(manifold.vknot).size()-1;
  std::printf("Default mesh has {%lu, %lu} elements\n",elu,elv);
  double A = 0.0;
  double Ael = 0.0;
  auto area = [](auto const &p) {return 1.0; };
  
  for (size_t j = 0; j < elv; j++)
  {
    for (size_t i = 0; i < elu; i++)
    {
      std::printf("\nIntegrating element(%lu,%lu)\n",i,j);
      std::cout << "updating element mesh..." << std::endl;
      mapper.updateElementMesh(i,j);
      std::cout << "updating integration points..." << std::endl;
      for (auto &p : ip) p.update();
      std::cout << "integrating weak forms..." << std::endl;
      Ael = 0.0;
      quadrature::gauss(ip,Ael,area);
      A += Ael;
    }
  }

  // double Aan = 3.0*M_PI;
  double Aan = 3.0*M_PI;
  double error = std::abs(Aan-A);
  std::printf("Computed area: %10.4g\nanalytc: %10.4g\nerror: %10.4g\n", A, Aan, error);
}

template<typename Solid>
void testGeometricBookKeeping(Solid const &solid)
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
  GeometricDofManager geom;
  geom.addShape(&solid);

  auto fake = nullptr;
  std::cout << geom.idsForShape(&solid) << std::endl;
  std::cout << geom.idsForShape(fake) << std::endl;

  // // get an iso surface
  NurbsSurface manifold; spline_ops::getIsoSurface(0.0,2,solid,manifold);
  geom.addShape(&manifold);
  std::cout << geom.idsForShape(&manifold) << std::endl;
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

int main(int argc, char **argv)
{
  NurbsSolid solid;
  elbow(1,2,3,solid);
  spline_ops::midpointRefinement(2,2,solid);
  // testGeometricBookKeeping(solid);
  // testSurfaceIntegration(solid);
  simulation(solid);

  return 0;
}
