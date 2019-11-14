#include "ManifoldElementMapper.h"

#include "splines/Nurbs.h"

#include "ParametricMesh.h"
#include "ShapeFunctions.h"


ManifoldElementMapper::
ManifoldElementMapper(NurbsSurface const &manifold)
  : ElementMapperBase()
  , _manifold(manifold)
  , _parametricMesh(iga::parametricMesh(_manifold))
{}

Eigen::RowVectorXd
ManifoldElementMapper::
_shape(double x1, double x2, double x3) const
{
  return iga::ShapeFunctions(x1,x2,_manifold);
}

Eigen::MatrixXd
ManifoldElementMapper::
_grad(double x1, double x2, double x3) const
{
  Eigen::MatrixXd const jacobian = _physicalJacobian(x1,x2,x3);
  Eigen::MatrixXd const dN = iga::ShapeFunctionDerivatives(x1,x2,_manifold);
  
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

  Eigen::MatrixXd    const Q = convert::to<Eigen::MatrixXd>(_manifold.Q);
  Eigen::MatrixXd    const jac = dN*Q;
  Eigen::RowVector3d const dxds   = jac.row(0);
  Eigen::RowVector3d const dxdt   = jac.row(1);

  Eigen::RowVector3d e1 = dxds.normalized();
  Eigen::RowVector3d e2 = dxdt.normalized();

  Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
  Eigen::Vector3d xc; xc << Q.col(0).mean(), Q.col(1).mean(), Q.col(2).mean();
  R.row(0) = e1;
  R.row(1) = e2;

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
  Eigen::MatrixXd const dN       = iga::ShapeFunctionDerivatives(x1,x2,_manifold);
  Eigen::MatrixXd const lQ       = computeLocalCoordinates(x1,x2,dN);
  Eigen::MatrixXd const jacobian = dN*lQ;

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
  _elementMesh = iga::parametricElementMesh(i,j,_manifold,_parametricMesh); 
}

