#include "splines/Algorithms.h"
#include "splines/BSpline.h"
#include "splines/Nurbs.h"
#include "splines/DegreeReduction.h"
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

void beam(NurbsSurface &surf, double width=1.0, double height=1.0)
{
  using namespace vector_ops;

  int p = 2;
  int q = 2;
  std::vector<double> uknot{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  algo::normalizeKnot(uknot);
  std::vector<double> vknot{0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  algo::normalizeKnot(vknot);
  typename NurbsCurve::matrix cpts{
    {0        , 0.0},
    {width/2.0, 0.0},
    {width    , 0.0},

    {0        , height/2.0},
    {width/2.0, height/2.0},
    {width    , height/2.0},

    {0        , height},
    {width/2.0, height},
    {width    , height}
  };

  surf.p = p;
  surf.q = q;
  surf.uknot = uknot;
  surf.vknot = vknot;
  surf.weights = std::vector<double>(cpts.size(),1.0);
  surf.Q = cpts;
}

void patch(NurbsSurface &surf, double radius=1.0, double width=2.0)
{
  using namespace vector_ops;
  double r1 = radius;

  int p = 2;
  int q = 1;
  std::vector<double> knot{0.0, 0.0, 0.0, 1, 1, 1};
  algo::normalizeKnot(knot);
  typename NurbsCurve::matrix cpts{
    {0.0, r1 },
    { r1, r1 },
    { r1 , 0.0},

    {0.0,    width},
    { width, width},
    { width, 0.0  },
  };

  static double const fact{1.0/std::sqrt(2.0)};
  decltype(knot) weights
  {
    1.0, fact, 1.0, 1.0, fact, 1.0
  };

  surf.p = p;
  surf.q = q;
  surf.uknot = knot;
  surf.vknot = {0,0,1,1};
  surf.weights = weights;
  surf.Q = cpts;
}

struct Plane
{
  std::vector<double> x;
  std::vector<double> n;
};

template<typename Shape>
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
  }

  return L;
}

template<typename Shape>
void simulation(Shape const &surf, Shape const &proj)
{
  std::cout << "\n+++ (" << __LINE__ << ") Enter: " << __PRETTY_FUNCTION__ << std::endl;
// #################################
// Operators 
// #################################
  auto Bop = [](auto const &p)->Eigen::MatrixXd
  {
    auto const grad = p.mapper.grad(p.para);
    auto const dof = grad.cols();
    Eigen::MatrixXd b(4,2*dof); b.setZero();

    for (int i = 0; i < dof; i++)
    {
      b(0,i*2)   = grad(0,i);
      b(1,i*2+1) = grad(1,i);

      b(3,i*2)   = grad(1,i);
      b(3,i*2+1) = grad(0,i);
    }
    return b;
  };

  auto Nop = [](auto const &p)->Eigen::MatrixXd
  {
    auto const shape = p.mapper.shape(p.para);
    auto const dof = shape.size();
    Eigen::MatrixXd N(2,2*dof); N.setZero();

    for (int i = 0; i < dof; i++)
    {
      N(0,i*2)   = shape[i];
      N(1,i*2+1) = shape[i];
    }
    return N;
  };

// #################################
// Material 
// #################################
  double const E = 1.0;
  double const nu = 0.4999999;
  // double const rho = 2700.0;
  auto matl = [E,nu](auto const &p)->Eigen::MatrixXd
  {
    Eigen::MatrixXd D(4,4); D.setZero();
    auto const mu = E/(2*(1+nu));
    auto const eye = Eigen::Matrix3d::Identity();
    double const ik[] = {0,1,2,0};
    double const jl[] = {0,1,2,1};
    for (size_t a = 0; a < 4; a++)
    {
      auto const i = ik[a];
      auto const j = jl[a];
      for (size_t b = a; b < 4; b++)
      {
        auto const k = ik[b];
        auto const l = jl[b];
        D(a,b) = mu*(eye(i,k)*eye(j,l) + eye(i,l)*eye(j,k));
        if (a != b) D(b,a) = D(a,b);
      }
    }
    return D;
  };

// #################################
// Weak forms 
// #################################
  // Vector2d accel; accel << 0.0, -9.81;
  // auto gravity = [&Nop, &rho, &accel](auto const &p)->Eigen::VectorXd
  // {
  //   auto const n = Nop(p);
  //   return rho*(n.transpose()*accel);
  // };

  auto stiffness = [&Bop,&matl](auto const &p)->Eigen::MatrixXd
  {
    auto const b = Bop(p);
    auto const d = matl(p);
    return b.transpose()*(d*b);
  };

  // set the up Simulation using the defined geometry
  GeometricDofManager geom;
  geom.addShape(&surf);

  // create mapper
  // // Compute stiffness
  size_t ndof = 2;
  size_t dof = ndof*geom.ctrlPoints.size();
  Eigen::MatrixXd K(dof,dof); K.setZero();
  Eigen::VectorXd f(dof); f.setZero();

  // isochoric stiffness
  {
    PlaneElementMapper mapper(surf);
    auto ip = iga::integrationPoints(mapper,surf.p,surf.q);

    // get assembly matrix for shape
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
        std::printf("\nIntegrating element(%lu,%lu)\n",i,j);
        std::cout << "updating element mesh..." << std::endl;
        mapper.updateElementMesh(i,j);
        std::cout << "updating integration points..." << std::endl;
        for (auto &p : ip) p.update();
        std::cout << "integrating weak forms..." << std::endl;
        Kel.setZero();
        quadrature::gauss(ip,Kel,stiffness);
        K += L*Kel*L.transpose();
      }
    }
  }

  // volumetric stiffness 
  {
    double const kb = E/(3*(1-2*nu));

    PlaneElementMapper rmapper(proj);

    auto rNop = [&rmapper](auto const &p)->Eigen::MatrixXd
    {
      auto const shape = rmapper.shape(p.para);
      return shape;
    };

    Eigen::VectorXd eye(4); eye << 1.0, 1.0, 1.0, 0.0;

    auto pBop = [&Bop, &rNop, eye](auto const &p)->Eigen::MatrixXd
    {
      auto const N = rNop(p);
      auto const b = Bop(p);
      return N.transpose()*(eye.transpose()*b);
    };

    auto Dvol = [&rNop,kb](auto const &p)->Eigen::MatrixXd
    {
      Eigen::MatrixXd const N = rNop(p);
      Eigen::MatrixXd const M = N.transpose()*N;
      Eigen::MatrixXd iM(Eigen::MatrixXd::Zero(M.rows(),M.cols()));
      auto sum = M.rowwise().sum();
      for (int i = 0; i < M.rows(); i++) {
        iM(i,i) = algo::equal(sum[i],0.0) ? 0.0 : 1.0/sum[i];
      }
      return kb*iM;
    };

    auto proj_stiffness = [&pBop,&Dvol](auto const &p)->Eigen::MatrixXd
    {
      auto const b = pBop(p);
      auto const d = Dvol(p);
      return b.transpose()*(d*b);
    };

    PlaneElementMapper mapper(surf);
    auto ip = iga::integrationPoints(mapper,proj.p,proj.q);

    // get assembly matrix for shape
    auto const L = assemblyMatrix(geom.idsForShape(&surf),ndof,geom.ctrlPoints.size());

    // integrate
    auto elu = iga::meshFromSpan(surf.uknot).size()-1;
    auto elv = iga::meshFromSpan(surf.vknot).size()-1;
    std::printf("\nDefault mesh has {%lu, %lu} elements\n",elu,elv);

    size_t edof = ndof*surf.Q.size();
    Eigen::MatrixXd Kel(edof,edof);

    for (size_t j = 0; j < elv; j++)
    {
      for (size_t i = 0; i < elu; i++)
      {
        std::printf("\nIntegrating element(%lu,%lu)\n",i,j);
        std::cout << "updating element mesh..." << std::endl;
        mapper.updateElementMesh(i,j);
        rmapper.updateElementMesh(i,j);
        std::cout << "updating integration points..." << std::endl;
        for (auto &p : ip) p.update();
        std::cout << "integrating weak forms..." << std::endl;
        Kel.setZero();
        quadrature::gauss(ip,Kel,proj_stiffness);
        K += L*Kel*L.transpose();
      }
    }
  }

// #################################
// Pressure Boundary condition
// #################################
  // define pressure form
  double papp = E/3.0;
  {
    auto pressure = [&Nop,papp](auto const &p)->Eigen::VectorXd
    {
      Eigen::MatrixXd const n = Nop(p);
      Eigen::VectorXd const normal = p.mapper.normal(p);

      return -papp*n.transpose()*normal;
    };

    // get the boundary fromt the surface
    NurbsCurve curve; spline_ops::getIsoCurve(0.0,1,surf,curve);
    geom.addShape(&curve);

    CurveElementMapper mapper(curve);
    // get some integration points for the manifold
    auto ip = iga::integrationPoints(mapper,curve.p);
    // get assembly matrix for shape
    auto const L = assemblyMatrix(geom.idsForShape(&curve),ndof,geom.ctrlPoints.size());
    // get the number of elements
    auto elu = iga::meshFromSpan(curve.knot).size()-1;
    std::printf("Default boundary mesh has {%lu} elements\n",elu);
    // integrate
    size_t edof = ndof*curve.Q.size();
    Eigen::VectorXd fel(edof);
    
    for (size_t i = 0; i < elu; i++)
    {
      std::printf("\nIntegrating boundary element(%lu)\n",i);
      std::cout << "updating boundary element mesh..." << std::endl;
      mapper.updateElementMesh(i);
      std::cout << "updating boundary integration points..." << std::endl;
      for (auto &p : ip) p.update();
      std::cout << "integrating boundary weak forms..." << std::endl;
      fel.setZero();
      quadrature::gauss(ip,fel,pressure);
      f += L*fel;
    }
    std::cout << "total load applied: " << f.sum() << std::endl;
  }

// #################################
// Dirichlet Boundary conditions
// #################################
  // compute constraints
  std::vector<size_t> cdof;
  // fix at the top
  Plane plane;
  plane.x = {0,0};
  plane.n = {-1,0};
  // x- plane
  auto const xsymm = pointsOnPlane<Shape>(geom,plane);
  for (auto const &idx : xsymm)
  {
    cdof.push_back(ndof*idx);
  }

  plane.n = {0,-1};
  auto const ysymm = pointsOnPlane<Shape>(geom,plane);
  for (auto const &idx : ysymm)
  {
    cdof.push_back(ndof*idx+1);
  }
  std::printf("Applying boundary conditions to %lu dof\n",cdof.size());

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
// solution 
// #################################
  auto AnalyticSolution = [E,nu,papp]()
  {
    auto Ri = 1.0;
    auto Ro = 4.0;
    auto A = ((1.0 + nu)*(Ri*Ro)*(Ri*Ro))/(E*(Ro*Ro - Ri*Ri));
    auto B = papp/Ri;
    auto C = (1.0-2.0*nu)*(papp*Ri)/(Ro*Ro);
    auto ur = A*(B + C);
    return ur;
  };
  std::cout << "**** Analytic solution = " << AnalyticSolution() << std::endl;

// #################################
// write solution to file and plot
// #################################
  // reshape solution
  auto const uv = Eigen::Map<const Eigen::MatrixXd>(u.data(),ndof,dof/ndof).transpose();
  // compute disp mag and deform the mesh
  auto const sids = geom.idsForShape(&surf);
  Eigen::VectorXd  umag(sids.size());
  Eigen::VectorXd vdisp(sids.size());

  double scale = 1.00;
  Shape deformed = surf;
  for (size_t i = 0; i < sids.size(); i++)
  {
    auto const x = convert::to<std::vector<double>>(Eigen::RowVectorXd(uv.row(sids[i])));
    umag[i]  = norm(x);
    vdisp[i] = x[1];
    deformed.Q[i][0] += scale*x[0];
    deformed.Q[i][1] += scale*x[1];
  }

  std::string file("output/nurbs_stress.txt");
  IO::writeSolutionToFile(deformed,umag,file,20,10);
  std::system(std::string(python + "python/plot_planar.py " + file).c_str());

  // check along v=w=const
  NurbsSurface wsolution; IO::geometryWithSolution(surf,umag,wsolution);
  size_t N = 100;
  std::vector<double> iso(N+1); std::iota(iso.begin(),iso.end(),0); iso/= N;
  double v = 0.0;
  std::vector<double> isoline;
  for (auto x : iso)
  {
    auto s = spline_ops::SurfacePoint(x,v,wsolution);
    isoline.push_back(s.back());
  }

  std::string isofile("output/nurbs_surface_iso.txt");
  IO::writeXYdata(iso,isoline,isofile);
  std::system(std::string(python + "python/plot_xy.py " + isofile).c_str());
  std::cout << "--- (" << __LINE__ << ") Exit: " << __PRETTY_FUNCTION__ << "\n" << std::endl;
}

template<typename Surface>
void SurfaceTest(Surface const &surf)
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  using namespace vector_ops;

  int N = 5;
  double u,v,du,dv;
  u = 0.0; v = 0.0; du = 1.0/N; dv = 1.0/N;
  std::vector<std::vector<double>> p, pu, pv;
  auto const Q = convert::to<Eigen::MatrixXd>(surf.Q);
  for (int i = 0; i <= N; i++)
  {
    u = 0.0;
    for (int i = 0; i <= N; i++)
    {
      p.push_back (spline_ops::SurfacePoint(u,v,surf));
      Eigen::MatrixXd const dN = iga::ShapeFunctionDerivatives(u,v,surf);
      Eigen::MatrixXd const grad = dN*Q;
      pu.push_back(convert::to<std::vector<double>>(Eigen::RowVectorXd(grad.row(0))));
      pv.push_back(convert::to<std::vector<double>>(Eigen::RowVectorXd(grad.row(1))));
      // pu.push_back(spline_ops::SurfaceDerivative(u,v,1,0,surf));
      // pv.push_back(spline_ops::SurfaceDerivative(u,v,1,1,surf));
      u += du; 
    }
    v += dv;
  }

  std::string file("output/nurbs_surface.txt");
  std::string uvec_file("output/nurbs_surface_du.txt");
  std::string vvec_file("output/nurbs_surface_dv.txt");
  spline_ops::writeVectorData(p,pu,uvec_file, true, 0.2);
  spline_ops::writeVectorData(p,pv,vvec_file, true, 0.2);
  spline_ops::writeToFile(surf,file,10,10);
  std::system(std::string(python + "python/plot_planar.py " + file + " " + uvec_file + " " + vvec_file).c_str());
}

template<typename Shape>
void reducedMesh(Shape const &shape)
{
  // SurfaceTest(shape);

  auto reduce = [&shape]()
  {
    Shape cpy;
    cpy.p = shape.p-1;
    cpy.q = shape.q-1;
    cpy.uknot.assign(shape.uknot.begin()+1,shape.uknot.end()-1);
    cpy.vknot.assign(shape.vknot.begin()+1,shape.vknot.end()-1);
    cpy.Q = shape.Q;
    cpy.weights = shape.weights;
    return cpy;
  };

  auto cpy = reduce();
  std::cout << shape << std::endl;
  std::cout << cpy << std::endl;
  // SurfaceTest(cpy);

}

template<typename Shape>
bool isReducible(Shape const &shape)
{
  std::cout << shape << std::endl;
  std::cout << "reduce in u..." << std::endl;
  Shape cpy1 = shape;
  spline_ops::insertKnot(0.5,2,0,cpy1);
  std::cout << cpy1 << std::endl;
  spline_ops::reduce(0,cpy1,1);
  std::cout << cpy1 << std::endl;
  std::cout << "\n\nreduce in v..." << std::endl;
  Shape cpy2 = shape;
  std::cout << cpy2 << std::endl;
  spline_ops::insertKnot(0.5,2,1,cpy2);
  std::cout << cpy2 << std::endl;
  spline_ops::reduce(1,cpy2,1);
  std::cout << cpy2 << std::endl;
  return false;
}

std::vector<std::pair<double,int>> knotMultiplicity(std::vector<double> const &knot)
{
  using pair_t = std::pair<double,int>;
  std::vector<pair_t> mult;
  if (knot.empty()) return mult;

  int m = 1;
  for (size_t i = 1; i < knot.size(); i++)
  {
    if (algo::equal(knot[i],knot[i-1])) m++;
    else
    {
      mult.push_back(std::make_pair(knot[i-1],m));
      m = 1;
    }
  }

  mult.push_back(std::make_pair(knot.back(),m));
  return mult;
}

template<typename Shape>
void makeCk(Shape &shape, int k, int direction)
{

  auto mult = knotMultiplicity((direction == 1) ? shape.vknot : shape.uknot);
  auto const order = (direction == 1) ?  shape.q : shape.p;
  for (auto const &p : mult)
  {
    auto diff = order - p.second;
    if (diff > k) spline_ops::insertKnot(p.first,diff-k,direction,shape);
  }
}

template<typename Shape>
void createMesh(Shape &shape, Shape &proj)
{
  std::cout << "\nInitial shape...." << std::endl;
  // std::cout << shape << std::endl;
  // insert some knots
  std::cout << "\n\np-refinement..." << std::endl;
  spline_ops::elevate(2,0,shape);
  spline_ops::elevate(2,1,shape);
  // std::cout << shape << std::endl;

  std::cout << "\n\nh-refinement..." << std::endl;
  spline_ops::midpointRefinement(1,0,shape);
  spline_ops::midpointRefinement(2,1,shape);
  // std::cout << shape << std::endl;

  std::cout << "\n\nmake ck..." << std::endl;
  int continuity = 0;
  makeCk(shape,continuity,0);
  makeCk(shape,continuity,1);
  // std::cout << shape << std::endl;

  std::cout << "\n\nprojection..." << std::endl;
  proj = shape;
  spline_ops::reduce(0,proj,1);
  spline_ops::reduce(1,proj,1);
  // std::cout << proj << std::endl;
  // 
  std::string file("output/foo.txt");
  spline_ops::writeToFile(shape,file,10,10);
  std::system(std::string(python + "python/plot_planar.py " + file).c_str());
  spline_ops::writeToFile(proj,file,10,10);
  std::system(std::string(python + "python/plot_planar.py " + file).c_str());
}

int main(int argc, char **argv)
{
  NurbsSurface surf,proj;
  patch(surf,1.0,4.0);
  createMesh(surf,proj);
  simulation(surf,proj);
  return 0;
}
