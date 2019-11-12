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
#include "iga/PlaneElementMapper.h"
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

void beam(NurbsSurface &surf, double width=1.0, double height=1.0)
{
  using namespace vector_ops;

  int p = 2;
  int q = 2;
  std::vector<double> uknot{0.0, 0.0, 0.0,1.0, 1.0, 1.0};
  algo::normalizeKnot(uknot);
  std::vector<double> vknot{0.0, 0.0, 0.0, 1.0,1.0, 1.0};
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
  double th = width-radius;
  double r1 = radius;
  double r2 = radius + th/2;

  int p = 2;
  int q = 2;
  std::vector<double> knot{0.0, 0.0, 0.0, 1, 1, 1};
  algo::normalizeKnot(knot);
  typename NurbsCurve::matrix cpts{
    {0.0,    width},
    {-width, width},
    {-width, 0.0  },
    
    {0.0, r2 },
    {-r2, r2 },
    {-r2, 0.0},

    {0.0, r1 },
    {-r1, r1 },
    {-r1, 0.0},
  };

  static double const fact{1.0/std::sqrt(2.0)};
  decltype(knot) weights
  {
    1.0, fact, 1.0, 1.0, fact ,1.0, 1.0, fact, 1.0
  };

  surf.p = p;
  surf.q = q;
  surf.uknot = knot;
  surf.vknot = knot;
  surf.weights = weights;
  surf.Q = cpts;
}

struct Plane
{
  std::vector<double> x;
  std::vector<double> n;
};


template<typename Shape>
class SimulationGeometry 
{
public:
  SimulationGeometry() = default;
  virtual ~SimulationGeometry() = default;

  bool addShape(Shape const *shape)
  {
    if (!shape) return false;

    if (std::find(shapes.begin(), shapes.end(),shape) == shapes.end())
    {
      shapes.push_back(shape);
      addNewCtrlPoints();
    } else {
      return false;
    }

    return true;
  }

  std::vector<size_t> const &idsForShape(Shape const *shape)
  {
    static std::vector<size_t> const zero;
    if (!shape) return zero;

    auto it = std::find(shapes.begin(),shapes.end(),shape);
    if (it != shapes.end())
      return sids[std::distance(shapes.begin(),it)];
    return zero;
  }

  std::vector<Shape const*> shapes; 
  std::vector<size_t> ids;
  std::vector<std::vector<size_t>> sids;
  std::vector<std::vector<double>> ctrlPoints;

private:
  void addNewCtrlPoints()
  {
    std::vector<size_t> tmp;
    for (auto const &x : shapes.back()->Q)
    {
      auto it = std::find(ctrlPoints.begin(), ctrlPoints.end(), x);
      if (it == ctrlPoints.end())
      {
        auto id = ctrlPoints.size();
        ids.push_back(id);
        tmp.push_back(id);
        ctrlPoints.push_back(x);
      } else {
        auto id = std::distance(ctrlPoints.begin(),it);
        tmp.push_back(id);
      }
    }

    sids.push_back(tmp);
  }
};

template<typename Shape>
std::vector<size_t> pointsOnPlane(SimulationGeometry<Shape> const &geom, Plane const &plane)
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
  double const E = 70e9;
  double const nu = 0.29;
  double const rho = 2700.0;
  auto matl = [E,nu](auto const &p)->Eigen::MatrixXd
  {
    Eigen::MatrixXd D(4,4); D.setZero();
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

    return D;
  };

// #################################
// Weak forms 
// #################################
  Vector2d accel; accel << 0.0, -9.81;
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

  // set the up Simulation using the defined geometry
  SimulationGeometry<Solid> geom;
  geom.addShape(&solid);

  // create mapper
  PlaneElementMapper mapper(solid);
  auto ip = iga::integrationPoints(mapper,solid.p,solid.q);

  // // Compute stiffness
  size_t ndof = 2;
  size_t dof = ndof*geom.ctrlPoints.size();
  size_t edof = ndof*solid.Q.size();
  Eigen::MatrixXd K(dof,dof); K.setZero();
  Eigen::VectorXd f(dof); f.setZero();

  Eigen::MatrixXd Kel(edof,edof);
  Eigen::VectorXd fel(edof);

  // get assembly matrix for shape
  auto const L = assemblyMatrix(geom.idsForShape(&solid),ndof,geom.ctrlPoints.size());

  // integrate 
  auto elu = iga::meshFromSpan(solid.uknot).size()-1;
  auto elv = iga::meshFromSpan(solid.vknot).size()-1;
  std::printf("Default mesh has {%lu, %lu} elements\n",elu,elv);
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
      Kel.setZero(); fel.setZero();
      quadrature::gauss(ip,Kel,stiffness);
      quadrature::gauss(ip,fel,gravity);
      K += L*Kel*L.transpose();
      f += L*fel;
    }
  }

// #################################
// Boundary conditions
// #################################
  // compute constraints
  // fix at the right
  Plane plane;
  plane.x = {0,0,0};
  plane.n = {1,0,0};
  // x- plane
  auto const fixed = pointsOnPlane<Solid>(geom,plane);
  std::printf("Applying boundary conditions for fixed surface to %lu cpts\n",fixed.size());

  size_t C = ndof*fixed.size();
  Eigen::MatrixXd Kc(C,dof); Kc.setZero();
  Eigen::VectorXd Rc(C); Rc.setZero();

  size_t k = 0;
  for (auto const &idx : fixed)
  {
    Kc(k++,ndof*idx)   = 1;
    Kc(k++,ndof*idx+1) = 1;
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
  auto const uv = Eigen::Map<const Eigen::MatrixXd>(u.data(),ndof,dof/ndof).transpose();
  // compute disp mag and deform the mesh
  auto const sids = geom.idsForShape(&solid);
  Eigen::VectorXd  umag(sids.size());
  Eigen::VectorXd vdisp(sids.size());

  double scale = 1.0;
  Solid deformed = solid;
  for (size_t i = 0; i < sids.size(); i++)
  {
    auto const x = convert::to<std::vector<double>>(Eigen::RowVectorXd(uv.row(sids[i])));
    umag[i]  = norm(x);
    vdisp[i] = x[1];
    deformed.Q[i][0] += scale*x[0];
    deformed.Q[i][1] += scale*x[1];
  }

  std::string file("output/nurbs_stress.txt");
  IO::writeSolutionToFile(deformed,vdisp,file,20,10);
  std::system(std::string(python + "python/plot_planar.py " + file).c_str());

  // check along v=w=const
  NurbsSurface wsolution; IO::geometryWithSolution(solid,vdisp,wsolution);
  size_t N = 100;
  std::vector<double> iso(N+1); std::iota(iso.begin(),iso.end(),0); iso/= N;
  double v = 1.0;
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

int main(int argc, char **argv)
{
  NurbsSurface surf;
  patch(surf,1.0,2.0);
  // beam(surf,4.0,2.0);
  SurfaceTest(surf);
  // std::string file("output/foo.txt");
  // spline_ops::writeToFile(surf,file,4,10);
  // std::system(std::string(python + "python/plot_planar.py " + file).c_str());
  // simulation(surf);
  return 0;
}
