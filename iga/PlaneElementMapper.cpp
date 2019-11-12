#include "PlaneElementMapper.h"

#include "splines/Nurbs.h"

#include "ParametricMesh.h"
#include "ShapeFunctions.h"


PlaneElementMapper::
PlaneElementMapper(NurbsSurface const &surface)
  : ElementMapperBase()
  , _surface(surface)
  , _parametricMesh(iga::parametricMesh(_surface))
{}

Eigen::RowVectorXd
PlaneElementMapper::
_shape(double x1, double x2, double x3) const
{
  return iga::ShapeFunctions(x1,x2,_surface);
}

Eigen::MatrixXd
PlaneElementMapper::
_grad(double x1, double x2, double x3) const
{
  Eigen::MatrixXd const jacobian = _physicalJacobian(x1,x2,x3);
  Eigen::MatrixXd const dN = iga::ShapeFunctionDerivatives(x1,x2,_surface);

  return jacobian.inverse()*dN;
}

Eigen::MatrixXd
PlaneElementMapper::
_Grad(double x1, double x2, double x3) const
{
  return _grad(x1,x2,x3);
}

Eigen::MatrixXd
PlaneElementMapper::
_parametricJacobian(double x1, double x2, double x3) const
{
  Eigen::MatrixXd const dN       = iga::parametricShapeFunctionDerivatives(x1,x2);
  Eigen::MatrixXd const jacobian = dN*_elementMesh;

  return jacobian.block(0,0,2,2);
}

Eigen::MatrixXd
PlaneElementMapper::
_physicalJacobian(double x1, double x2, double x3) const
{
  Eigen::MatrixXd const Q        = convert::to<Eigen::MatrixXd>(_surface.Q);
  Eigen::MatrixXd const dN       = iga::ShapeFunctionDerivatives(x1,x2,_surface);
  Eigen::MatrixXd const jacobian = dN*Q;
  return jacobian.block(0,0,2,2);
}

void
PlaneElementMapper::
_mapIntegrationPoint(double  x1, double  x2, double  x3,
                     double &p1, double &p2, double &p3) const

{
  Eigen::RowVectorXd const N = iga::parametricShapeFunction(x1,x2);
  Eigen::RowVectorXd const mapped{N*_elementMesh};

  p1 = mapped[0];
  p2 = mapped[1];
}

void
PlaneElementMapper::
_updateElementMesh(size_t i, size_t j, size_t k)
{
  _elementMesh = iga::parametricElementMesh(i,j,_surface,_parametricMesh); 
}

