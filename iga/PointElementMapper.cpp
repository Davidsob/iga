#include "PointElementMapper.h"

#include "splines/GeometricObject.h"
#include "splines/BSpline.h"
#include "splines/utils/VectorOperations.h"

#include "GeometricDofManager.h"
#include "ParametricMesh.h"
#include "ShapeFunctions.h"


PointElementMapper::
PointElementMapper(BSplinePoint const &point)
  : ElementMapperBase()
  , _point(point)
  , _parametricMesh(Eigen::MatrixXd::Zero(1,1))
{
  auto &mgr = GeometricDofManager::instance();
  mgr.addShape(&_point);
  _dof = mgr.idsForShape(&_point);
}

GeometricObject const *
PointElementMapper::geometry() const
{ 
  return dynamic_cast<GeometricObject const *>(&_point);
}

Eigen::RowVectorXd
PointElementMapper::
_shape(double x1, double x2, double x3) const
{
  static Eigen::RowVectorXd one(1); one << 1;
  return one;
}

Eigen::MatrixXd
PointElementMapper::
_grad(double x1, double x2, double x3) const
{
  static Eigen::MatrixXd const zero(Eigen::MatrixXd::Zero(1,1));
  return zero;
}

Eigen::MatrixXd
PointElementMapper::
_Grad(double x1, double x2, double x3) const
{
  return _grad(x1,x2,x3);
}

Eigen::MatrixXd
PointElementMapper::
_parametricJacobian(double x1, double x2, double x3) const
{
  static Eigen::MatrixXd one(1,1); one << 1;
  return one;
}

Eigen::MatrixXd
PointElementMapper::
_physicalJacobian(double x1, double x2, double x3) const
{
  static Eigen::MatrixXd one(1,1); one << 1;
  return one;
}

void
PointElementMapper::
_mapIntegrationPoint(double  x1, double  x2, double  x3,
                     double &p1, double &p2, double &p3) const

{
}

void
PointElementMapper::
_updateElementMesh(size_t i, size_t j, size_t k)
{
  auto &mgr = GeometricDofManager::instance();

  _elementMesh = _parametricMesh; 
  _dof = mgr.idsForShape(&_point);
  _coordinates = convert::to<Eigen::MatrixXd>(subVector(mgr.ctrlPoints,_dof));
}

