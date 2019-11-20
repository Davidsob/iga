#include "CurveElementMapper.h"

#include "splines/Nurbs.h"

#include "ParametricMesh.h"
#include "ShapeFunctions.h"


CurveElementMapper::
CurveElementMapper(NurbsCurve const &curve)
  : ElementMapperBase()
  , _curve(curve)
  , _parametricMesh(iga::parametricMesh(_curve))
{}

Eigen::RowVectorXd
CurveElementMapper::
_shape(double x1, double x2, double x3) const
{
  return iga::ShapeFunctions(x1,_curve);
}

Eigen::MatrixXd
CurveElementMapper::
_grad(double x1, double x2, double x3) const
{
  Eigen::MatrixXd const jacobian = _physicalJacobian(x1,x2,x3);
  Eigen::MatrixXd const dN = iga::ShapeFunctionDerivatives(x1,_curve);
  
  return jacobian.inverse()*dN;
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

// Eigen::MatrixXd
// CurveElementMapper::
// computeLocalCoordinates(double x1, Eigen::MatrixXd const &dN) const
// {
//   static Eigen::RowVector3d e3({0,0,1});

//   Eigen::MatrixXd    const Q = convert::to<Eigen::MatrixXd>(_curve.Q);
//   Eigen::MatrixXd    const jac = dN*Q;

//   auto const dim = Q.cols();
//   Eigen::MatrixXd R = Eigen::MatrixXd::Zero(dim,dim);

//   Eigen::Vector3d e1; e1.head(2) = jac.row(0).normalized(); e1[2] = 0.0;
//   Eigen::Vector3d e2 = e3.cross(e1); e2 = e2.normalized();

//   R.row(0) = e1.head(2);
//   // R.row(1) = e2.head(2);

//   Eigen::MatrixXd lQ(Q.rows(),1);
//   Eigen::VectorXd xc = _shape(x1,0,0)*Q;
//   for (int i = 0; i < Q.rows(); i++)
//   {
//     auto xg = Eigen::VectorXd(Q.row(i));
//     auto xl = R*(xg-xc)+xc;
//     lQ(i,0) = xl(0,0);
//   }
//   return lQ;
// }

Eigen::MatrixXd
CurveElementMapper::
computeLocalCoordinates(double x1, Eigen::MatrixXd const &dN) const
{
  Eigen::MatrixXd  const Q = convert::to<Eigen::MatrixXd>(_curve.Q);
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
  Eigen::MatrixXd const dN       = iga::ShapeFunctionDerivatives(x1,_curve);
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
  _elementMesh = iga::parametricElementMesh(i,_curve,_parametricMesh); 
}

Eigen::VectorXd
CurveElementMapper::
_normal(double x1) const
{
  static const double eps = 1e-6;
  static const Eigen::Vector3d e3{0,0,1};

  Eigen::MatrixXd const Q  = convert::to<Eigen::MatrixXd>(_curve.Q);
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
