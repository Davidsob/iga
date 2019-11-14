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

  NurbsSurface const &_manifold;
  Eigen::MatrixXd _parametricMesh;
  mutable Eigen::MatrixXd _elementMesh;
};

inline std::ostream &operator<<(std::ostream &os, ManifoldElementMapper const &mapper)
{
  os << "\nManifold Element Mapper:" << std::endl;
  os << "*** parametric mesh:\n" << mapper.parametricMesh() << std::endl;
  os << "*** element mesh:\n" << mapper.elementMesh() << std::endl;
  return os;
}
