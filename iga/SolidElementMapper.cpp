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
  Eigen::MatrixXd const jacobian = _physicalJacobian(x1,x2,x3);
  Eigen::MatrixXd const dN = iga::ShapeFunctionDerivatives(x1,x2,x3,_solid);
  
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
  Eigen::MatrixXd const Q        = convert::to<Eigen::MatrixXd>(_solid.Q);
  Eigen::MatrixXd const dN       = iga::ShapeFunctionDerivatives(x1,x2,x3,_solid);
  Eigen::MatrixXd const jacobian = dN*Q;

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
  _elementMesh = iga::parametricElementMesh(i,j,k,_solid,_parametricMesh); 
}

