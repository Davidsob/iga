#include "SolidElementMapper.h"

#include "splines/GeometricObject.h"
#include "splines/Nurbs.h"
#include "splines/utils/VectorOperations.h"

#include "GeometricDofManager.h"
#include "ParametricMesh.h"
#include "ShapeFunctions.h"


SolidElementMapper::
SolidElementMapper(NurbsSolid const &solid)
  : ElementMapperBase()
  , _solid(solid)
  , _parametricMesh(iga::parametricMesh(_solid))
{
  auto &mgr = GeometricDofManager::instance();
  mgr.addShape(&_solid);
  _dof = mgr.idsForShape(&_solid);
}

GeometricObject const *
SolidElementMapper::geometry() const
{ 
  return dynamic_cast<GeometricObject const *>(&_solid);
}

Eigen::RowVectorXd
SolidElementMapper::
_shape(double x1, double x2, double x3) const
{
  return iga::CompactShapeFunctions(x1,x2,x3,_solid);
}

Eigen::MatrixXd
SolidElementMapper::
_grad(double x1, double x2, double x3) const
{
  Eigen::MatrixXd const jacobian = _physicalJacobian(x1,x2,x3);
  Eigen::MatrixXd const dN = iga::CompactShapeFunctionDerivatives(x1,x2,x3,_solid);
  
  return jacobian.inverse()*dN;
}

Eigen::MatrixXd
SolidElementMapper::
_Grad(double x1, double x2, double x3) const
{
  return _grad(x1,x2,x3);
}

Eigen::MatrixXd
SolidElementMapper::
_parametricJacobian(double x1, double x2, double x3) const
{
  Eigen::MatrixXd const dN       = iga::parametricShapeFunctionDerivatives(x1,x2,x3);
  Eigen::MatrixXd const jacobian = dN*_elementMesh;

  return jacobian;
}

Eigen::MatrixXd
SolidElementMapper::
_physicalJacobian(double x1, double x2, double x3) const
{
  Eigen::MatrixXd const dN       = iga::CompactShapeFunctionDerivatives(x1,x2,x3,_solid);
  Eigen::MatrixXd const jacobian = dN*_coordinates;

  return jacobian;
}

void
SolidElementMapper::
_mapIntegrationPoint(double  x1, double  x2, double  x3,
                     double &p1, double &p2, double &p3) const

{
  Eigen::RowVectorXd const N = iga::parametricShapeFunction(x1,x2,x3);
  Eigen::RowVectorXd const mapped{N*_elementMesh};

  p1 = mapped[0];
  p2 = mapped[1];
  p3 = mapped[2];
}

void
SolidElementMapper::
_updateElementMesh(size_t i, size_t j, size_t k)
{
  auto &mgr = GeometricDofManager::instance();

  _elementMesh = iga::parametricElementMesh(i,j,k,_solid,_parametricMesh); 
  auto const uvw{_elementMesh.row(0)};
  auto const active{iga::ActiveControlPoints(uvw(0),uvw(1),uvw(2),_solid)};
  auto const sids{mgr.idsForShape(&_solid)};
  _dof = subVector(sids,active);
  _coordinates = convert::to<Eigen::MatrixXd>(subVector(mgr.ctrlPoints,_dof));
}

