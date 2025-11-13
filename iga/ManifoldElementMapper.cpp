#include "ManifoldElementMapper.h"

#include "splines/GeometricObject.h"
#include "splines/Nurbs.h"
#include "splines/utils/VectorOperations.h"

#include "GeometricDofManager.h"
#include "ParametricMesh.h"
#include "ShapeFunctions.h"


ManifoldElementMapper::
ManifoldElementMapper(NurbsSurface const &manifold)
  : ElementMapperBase()
  , _manifold(manifold)
  , _parametricMesh(iga::parametricMesh(_manifold))
{
  auto &mgr = GeometricDofManager::instance();
  mgr.addShape(&_manifold);
}

GeometricObject const *
ManifoldElementMapper::geometry() const
{ 
  return dynamic_cast<GeometricObject const *>(&_manifold);
}

Eigen::RowVectorXd
ManifoldElementMapper::
_shape(double x1, double x2, double x3) const
{
  return iga::CompactShapeFunctions(x1,x2,_manifold);
}

Eigen::MatrixXd
ManifoldElementMapper::
_grad(double x1, double x2, double x3) const
{
  Eigen::MatrixXd const jacobian = _physicalJacobian(x1,x2,x3);
  Eigen::MatrixXd const dN = iga::CompactShapeFunctionDerivatives(x1,x2,_manifold);
  
  return jacobian.inverse()*dN;
}

Eigen::MatrixXd
ManifoldElementMapper::
_Grad(double x1, double x2, double x3) const
{
  return _grad(x1,x2,x3);
}

Eigen::MatrixXd
ManifoldElementMapper::
_parametricJacobian(double x1, double x2, double x3) const
{
  Eigen::MatrixXd const dN       = iga::parametricShapeFunctionDerivatives(x1,x2);
  Eigen::MatrixXd const jacobian = dN*_elementMesh;

  return jacobian;
}

Eigen::MatrixXd
ManifoldElementMapper::
computeLocalCoordinates(double x1, double x2, Eigen::MatrixXd const &dN) const
{
  static const std::vector<double> n0{0,0,1};

  Eigen::MatrixXd    const &Q  = _coordinates;
  Eigen::MatrixXd    const jac = dN*Q;
  Eigen::RowVector3d const dxds   = jac.row(0);
  Eigen::RowVector3d const dxdt   = jac.row(1);

  Eigen::RowVector3d e1 = dxds.normalized();
  Eigen::RowVector3d e2 = dxdt.normalized();

  Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
  R.row(0) = e1;
  R.row(1) = e2;
  
  Eigen::Vector3d xc = Q.colwise().mean();

  Eigen::MatrixXd lQ(Q.rows(),2);
  for (int i = 0; i < Q.rows(); i++)
  {
    auto xg = Eigen::Vector3d(Q.row(i));
    auto xl = R*(xg-xc) + xc;
    lQ.row(i) = xl.head(2);
  }

  return lQ;
}

Eigen::MatrixXd
ManifoldElementMapper::
_physicalJacobian(double x1, double x2, double x3) const
{
  Eigen::MatrixXd const dN       = iga::CompactShapeFunctionDerivatives(x1,x2,_manifold);
  Eigen::MatrixXd const lQ       = computeLocalCoordinates(x1,x2,dN);
  Eigen::MatrixXd const jacobian = dN*lQ;
  // return covariantMetricTensor(std::vector<double>{x1,x2,x3}).block(0,0,2,2);
  return jacobian;
}

void
ManifoldElementMapper::
_mapIntegrationPoint(double  x1, double  x2, double  x3,
                     double &p1, double &p2, double &p3) const

{
  Eigen::RowVectorXd const N = iga::parametricShapeFunction(x1,x2);
  Eigen::RowVectorXd const mapped{N*_elementMesh};

  p1 = mapped[0];
  p2 = mapped[1];
}

void
ManifoldElementMapper::
_updateElementMesh(size_t i, size_t j, size_t k)
{
  auto &mgr = GeometricDofManager::instance();

  _elementMesh = iga::parametricElementMesh(i,j,_manifold,_parametricMesh); 
  auto const uv{_elementMesh.row(0)};
  auto const active{iga::ActiveControlPoints(uv(0),uv(1),_manifold)};
  auto const sids{mgr.idsForShape(&_manifold)};

  _dof = subVector(sids,active);
  _coordinates = convert::to<Eigen::MatrixXd>(subVector(mgr.ctrlPoints,_dof));
}

Eigen::VectorXd
ManifoldElementMapper::
_normal(double x1, double x2) const
{
  Eigen::MatrixXd const tangents = _tangents(x1,x2);
  Eigen::Vector3d const &e1 = tangents.row(0);
  Eigen::Vector3d const &e2 = tangents.row(1);
  Eigen::Vector3d const normal = e1.cross(e2);
  return normal.normalized();
}

Eigen::MatrixXd
ManifoldElementMapper::
_tangents(double x1, double x2) const
{
  Eigen::MatrixXd const &Q = _coordinates;
  Eigen::MatrixXd const  g = _grad(x1,x2,0);

  Eigen::MatrixXd tangents = g*Q;
  tangents.row(0) = tangents.row(0).normalized();
  tangents.row(1) = tangents.row(1).normalized();

  return tangents;
}

Eigen::Matrix3d
ManifoldElementMapper::
_localTransformation(double x1, double x2) const
{
  Eigen::MatrixXd const tangents = _tangents(x1,x2);
  Eigen::Vector3d const &e1 = tangents.row(0);
  Eigen::Vector3d const &e2 = tangents.row(1);
  Eigen::Vector3d const normal = e1.cross(e2);

  Eigen::Matrix3d R;
  R.block(0,0,2,3) = tangents;
  R.row(2) = normal;

  return R;
}

Eigen::Matrix3d
ManifoldElementMapper::
_covariantBasis(double x1, double x2) const
{
  Eigen::MatrixXd const &Q    = _coordinates;
  Eigen::MatrixXd const sgrad = iga::CompactShapeFunctionDerivatives(x1,x2,_manifold);

  Eigen::Matrix3d A;
  A.block(0,0,2,3) = sgrad*Q;

  Eigen::Vector3d N = A.row(0).cross(A.row(1));
  A.row(2) = N.normalized();

  return A.transpose();
}

