#include "SolidElementMapper.h"

#include "splines/Nurbs.h"

#include "ParametricMesh.h"
#include "ShapeFunctions.h"


SolidElementMapper::
SolidElementMapper(NurbsSolid const &solid)
  : ElementMapperBase()
  , _solid(solid)
  , _parametricMesh(iga::parametricMesh(_solid))
{}

Eigen::RowVectorXd
SolidElementMapper::
_shape(double x1, double x2, double x3) const
{
  return iga::ShapeFunctions(x1,x2,x3,_solid);
}

Eigen::MatrixXd
SolidElementMapper::
_grad(double x1, double x2, double x3) const
{
  auto const jacobian = _physicalJacobian(x1,x2,x3);
  auto const dN = iga::ShapeFunctionDerivatives(x1,x2,x3,_solid);
  return jacobian.inverse()*dN;
}

Eigen::MatrixXd
SolidElementMapper::
_Grad(double x1, double x2, double x3) const
{
  return _grad(x1,x2,x3);
}

Eigen::Matrix3d
SolidElementMapper::
_parametricJacobian(double x1, double x2, double x3) const
{
  auto const dN       = iga::parametricShapeFunctionDerivatives(x1,x2,x3);
  auto const jacobian = dN*_elementMesh;

  return jacobian.transpose();
}

Eigen::Matrix3d
SolidElementMapper::
_physicalJacobian(double x1, double x2, double x3) const
{
  auto const dN       = iga::ShapeFunctionDerivatives(x1,x2,x3,_solid);
  auto const jacobian = dN*convert::to<MatrixXd>(_solid.Q);

  return jacobian.transpose();
}

void
SolidElementMapper::
_mapIntegrationPoint(double  x1, double  x2, double  x3,
                     double &p1, double &p2, double &p3) const

{
  auto const N = iga::parametricShapeFunction(x1,x2,x3);
  Eigen::RowVectorXd const mapped{N*_elementMesh};

  p1 = mapped[0];
  p2 = mapped[1];
  p3 = mapped[2];
}

void
SolidElementMapper::
_updateElementMesh(size_t i, size_t j, size_t k)
{
  _elementMesh = iga::parametricElementMesh(i,j,k,_solid,_parametricMesh); 
}

