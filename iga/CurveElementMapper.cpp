#include "CurveElementMapper.h"

#include "splines/GeometricObject.h"
#include "splines/Nurbs.h"
#include "splines/utils/VectorOperations.h"

#include "GeometricDofManager.h"
#include "ParametricMesh.h"
#include "ShapeFunctions.h"


CurveElementMapper::
CurveElementMapper(NurbsCurve const &curve)
  : ElementMapperBase()
  , _curve(curve)
  , _parametricMesh(iga::parametricMesh(_curve))
{
  auto &mgr = GeometricDofManager::instance();
  mgr.addShape(&_curve);
}

GeometricObject const *
CurveElementMapper::geometry() const
{ 
  return dynamic_cast<GeometricObject const *>(&_curve);
}

Eigen::RowVectorXd
CurveElementMapper::
_shape(double x1, double x2, double x3) const
{
  return iga::CompactShapeFunctions(x1,_curve);
}

Eigen::MatrixXd
CurveElementMapper::
_grad(double x1, double x2, double x3) const
{
  Eigen::MatrixXd const jacobian = _physicalJacobian(x1,x2,x3);
  Eigen::MatrixXd const dN = iga::CompactShapeFunctionDerivatives(x1,_curve);
  
  return jacobian.inverse()*dN;
}

Eigen::MatrixXd
CurveElementMapper::
_grad2(double x1) const
{
  auto const Q = subVector(_curve.Q, _dof);
  Eigen::MatrixXd const jacobian = _physicalJacobian(x1,0,0);
  Eigen::MatrixXd dN2(1,Q.size());
  dN2.row(0) = iga::CompactShapeFunctionDerivative2(x1,_curve);

  auto ijac = jacobian.inverse();
  return ijac*ijac*dN2;
}

Eigen::MatrixXd
CurveElementMapper::
_Grad(double x1, double x2, double x3) const
{
  return _grad(x1,x2,x3);
}

Eigen::MatrixXd
CurveElementMapper::
_parametricJacobian(double x1, double x2, double x3) const
{
  Eigen::MatrixXd const dN       = iga::parametricShapeFunctionDerivatives(x1);
  Eigen::MatrixXd const jacobian = dN*_elementMesh;

  return jacobian;
}

Eigen::MatrixXd
CurveElementMapper::
computeLocalCoordinates(double x1, Eigen::MatrixXd const &dN) const
{
  Eigen::MatrixXd  const &Q = _coordinates;
  Eigen::MatrixXd  const jac = dN*Q;
  Eigen::VectorXd  const x = _shape(x1,0,0)*Q;

  Eigen::VectorXd e1 = jac.row(0).normalized();
  Eigen::MatrixXd lQ(Q.rows(),1);
  for (int i = 0; i < Q.rows(); i++)
  {
    auto xg = Eigen::VectorXd(Q.row(i));
    auto xl = (xg-x).dot(e1);
    lQ(i,0) = xl;
  }
  return lQ;
}

Eigen::MatrixXd
CurveElementMapper::
_physicalJacobian(double x1, double x2, double x3) const
{
  Eigen::MatrixXd const dN       = iga::CompactShapeFunctionDerivatives(x1,_curve);
  Eigen::MatrixXd const lQ       = computeLocalCoordinates(x1,dN);
  Eigen::MatrixXd const jacobian = dN*lQ;
  return jacobian;
}

void
CurveElementMapper::
_mapIntegrationPoint(double  x1, double  x2, double  x3,
                     double &p1, double &p2, double &p3) const

{
  Eigen::RowVectorXd const N = iga::parametricShapeFunction(x1);
  Eigen::RowVectorXd const mapped{N*_elementMesh};

  p1 = mapped[0];
}

void
CurveElementMapper::
_updateElementMesh(size_t i, size_t j, size_t k)
{
  auto &mgr = GeometricDofManager::instance();
  
  _elementMesh = iga::parametricElementMesh(i,_curve,_parametricMesh); 
  auto const u{_elementMesh.row(0)};
  auto const active{iga::ActiveControlPoints(u(0),_curve)};
  auto const sids{mgr.idsForShape(&_curve)};
  _dof = subVector(sids,active);
  _coordinates = convert::to<Eigen::MatrixXd>(subVector(mgr.ctrlPoints,_dof));
}

double
CurveElementMapper::
_curvature(double x1) const
{
  auto dt = spline_ops::CurveDerivatives(x1,2,_curve);
  auto num = norm(cross(dt[1],dt[2]));
  auto den = std::pow(norm(dt[1]),3);

  return algo::equal(den,0.0) ? 0.0 : num/den;
}

Eigen::VectorXd
CurveElementMapper::
_tangent(double x1) const
{
  Eigen::VectorXd const t   = Eigen::MatrixXd(_grad(x1,0,0)*_coordinates).row(0).normalized();

  return t.normalized();
}

Eigen::VectorXd
CurveElementMapper::
_normal(double x1) const
{
  static const double eps = 1e-6;
  static const Eigen::Vector3d e3{0,0,1};

  Eigen::MatrixXd const &Q  = _coordinates;
  Eigen::VectorXd const t   = Eigen::MatrixXd(_grad(x1,0,0)*Q).row(0).normalized();
  Eigen::VectorXd const tf  = Eigen::MatrixXd(_grad(x1+eps,0,0)*Q).row(0).normalized();

  Eigen::VectorXd N = (tf-t)/eps;
  if (algo::equal(N.norm(),0.0))
  {
    auto const dim = Q.cols();
    Eigen::Vector3d e1; e1.head(dim) = t;
    return e1.cross(e3).head(dim);
  } else {
    return N.normalized();
  }
}

Eigen::MatrixXd
CurveElementMapper::
_localTransformation(double const x1) const
{
  static const Eigen::Vector3d e3{0,0,1};

  Eigen::MatrixXd const &Q  = _coordinates;
  Eigen::VectorXd const t   = Eigen::MatrixXd(_grad(x1,0,0)*Q).row(0).normalized();

  auto const dim = Q.cols();
  Eigen::Vector3d e1; e1.head(dim) = t; 
  Eigen::RowVector3d e2 = e3.cross(e1);
  Eigen::MatrixXd R(dim,dim);
  R.row(0) = e1.head(dim);
  R.row(1) = e2.head(dim);

  return R;
}
