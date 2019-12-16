#pragma once

#include "ElementMapperBase.h"

class NurbsSurface;

class ManifoldElementMapper
 : public ElementMapperBase
{
public:

  explicit ManifoldElementMapper(NurbsSurface const &manifold);

  ~ManifoldElementMapper() = default;

  friend std::ostream &operator<<(std::ostream &os, ManifoldElementMapper const &mapper);

  Eigen::MatrixXd const &parametricMesh() const { return _parametricMesh; }
  Eigen::MatrixXd const &elementMesh()    const { return _elementMesh; }

  std::vector<size_t>const &dof() const override { return _dof;}

  GeometricObject const * geometry() const override;

  Eigen::MatrixXd const &coordinates() const override { return _coordinates; } 

  NurbsSurface const &manifold() const { return _manifold; }
  
  template<typename ParametricPoint>
  Eigen::VectorXd normal(ParametricPoint const &p) const
  {
    return _normal(p[0],p[1]);
  }

  template<typename ParametricPoint>
  Eigen::MatrixXd tangents(ParametricPoint const &p) const
  {
    return _tangents(p[0],p[1]);
  }

  template<typename ParametricPoint>
  Eigen::Matrix3d localTransformation(ParametricPoint const &p) const
  {
    return _localTransformation(p[0],p[1]);
  }

  template<typename ParametricPoint>
  Eigen::Matrix3d covariantBasis(ParametricPoint const &p) const
  {
    return _covariantBasis(p[0],p[1]);
  }

  template<typename ParametricPoint>
  Eigen::Matrix3d contravariantBasis(ParametricPoint const &p) const
  {
    return contravariantMetricTensor(p)*covariantBasis(p);
  }

  template<typename ParametricPoint>
  Eigen::Matrix3d covariantMetricTensor(ParametricPoint const &p) const
  {
    auto const A = covariantBasis(p);
    return A.transpose()*A;
  }

  template<typename ParametricPoint>
  Eigen::Matrix3d contravariantMetricTensor(ParametricPoint const &p) const
  {
    return covariantMetricTensor(p).inverse();
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

  Eigen::MatrixXd computeLocalCoordinates(double x1, double x2, Eigen::MatrixXd const &dN) const;
  Eigen::VectorXd _normal(double x1, double x2) const;
  Eigen::MatrixXd _tangents(double x1, double x2) const;
  Eigen::Matrix3d _localTransformation(double x1, double x2) const;
  Eigen::Matrix3d _covariantBasis(double x1, double x2) const;

  NurbsSurface const &_manifold;
  std::vector<size_t> _dof;
  Eigen::MatrixXd _parametricMesh;
  mutable Eigen::MatrixXd _elementMesh;
  mutable Eigen::MatrixXd _coordinates;
};

inline std::ostream &operator<<(std::ostream &os, ManifoldElementMapper const &mapper)
{
  os << "\nManifold Element Mapper:" << std::endl;
  os << "*** parametric mesh:\n" << mapper.parametricMesh() << std::endl;
  os << "*** element mesh:\n" << mapper.elementMesh() << std::endl;
  return os;
}
