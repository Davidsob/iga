#pragma once

#include "ElementMapperBase.h"
#include <iostream>

class NurbsSurface;

class PlaneElementMapper
 : public ElementMapperBase
{
public:

  explicit PlaneElementMapper(NurbsSurface const &surface);

  ~PlaneElementMapper() = default;

  friend std::ostream &operator<<(std::ostream &os, PlaneElementMapper const &mapper);

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
  NurbsSurface const &_surface;
  Eigen::MatrixXd _parametricMesh;
  mutable Eigen::MatrixXd _elementMesh;
};

inline std::ostream &operator<<(std::ostream &os, PlaneElementMapper const &mapper)
{
  os << "\nPlane Element Mapper:" << std::endl;
  os << "*** parametric mesh:\n" << mapper.parametricMesh() << std::endl;
  os << "*** element mesh:\n" << mapper.elementMesh() << std::endl;
  return os;
}
