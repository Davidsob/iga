#pragma once

#include "ElementMapperBase.h"

class NurbsCurve;

class CurveElementMapper
 : public ElementMapperBase
{
public:

  explicit CurveElementMapper(NurbsCurve const &curve);

  ~CurveElementMapper() = default;

  friend std::ostream &operator<<(std::ostream &os, CurveElementMapper const &mapper);

  Eigen::MatrixXd const &parametricMesh() const { return _parametricMesh; }
  Eigen::MatrixXd const &elementMesh()    const { return _elementMesh; }

  std::vector<size_t>const &dof() const override { return _dof;}

  GeometricObject const * geometry() const override;

  template<typename ParametricPoint>
  Eigen::VectorXd tangent(ParametricPoint const &p) const
  {
    return _tangent(p[0]);
  }

  template<typename ParametricPoint>
  Eigen::VectorXd normal(ParametricPoint const &p) const
  {
    return _normal(p[0]);
  }

  template<typename ParametricPoint>
  double curvature(ParametricPoint const &p) const
  {
    return _curvature(p[0]);
  }

  template<typename ParametricPoint>
  Eigen::MatrixXd grad2(ParametricPoint const &p) const
  {
    return _grad2(p[0]);
  }

  template<typename ParametricPoint>
  Eigen::MatrixXd localTransformation(ParametricPoint const &p) const
  {
    return _localTransformation(p[0]);
  }

protected:
  Eigen::RowVectorXd _shape(double x1, double x2, double x3)                  const override;
  Eigen::MatrixXd    _grad(double x1, double x2, double x3)                   const override;
  Eigen::MatrixXd    _Grad(double x1, double x2, double x3)                   const override;
  Eigen::MatrixXd    _parametricJacobian(double x1, double x2, double x3)     const override;
  Eigen::MatrixXd    _physicalJacobian(double x1, double x2, double x3)       const override;
  void               _mapIntegrationPoint(double  x1, double  x2, double  x3,
                                          double &p1, double &p2, double &p3) const override;
  void               _updateElementMesh(size_t i, size_t j, size_t k)               override;

private:

  Eigen::MatrixXd computeLocalCoordinates(double x1, Eigen::MatrixXd const &dN) const;
  Eigen::VectorXd _tangent(double x1) const;
  Eigen::VectorXd _normal(double x1) const;
  double          _curvature(double x1) const;
  Eigen::MatrixXd _grad2(double x1) const;
  Eigen::MatrixXd _localTransformation(double x1) const;

  NurbsCurve const &_curve;
  std::vector<size_t> _dof;
  Eigen::MatrixXd _parametricMesh;
  mutable Eigen::MatrixXd _elementMesh;
};

inline std::ostream &operator<<(std::ostream &os, CurveElementMapper const &mapper)
{
  os << "\nCurve Element Mapper:" << std::endl;
  os << "*** parametric mesh:\n" << mapper.parametricMesh() << std::endl;
  os << "*** element mesh:\n" << mapper.elementMesh() << std::endl;
  return os;
}
