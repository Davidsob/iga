#pragma once

#include "ElementMapperBase.h"

class BSplinePoint;

class PointElementMapper
 : public ElementMapperBase
{
public:

  explicit PointElementMapper(BSplinePoint const &point);

  ~PointElementMapper() = default;

  friend std::ostream &operator<<(std::ostream &os, PointElementMapper const &mapper);

  Eigen::MatrixXd const &parametricMesh() const { return _parametricMesh; }
  Eigen::MatrixXd const &elementMesh()    const { return _elementMesh; }

  std::vector<size_t>const &dof() const override { return _dof;}

  GeometricObject const * geometry() const override;

  Eigen::MatrixXd const &coordinates() const override { return _coordinates; } 

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
  BSplinePoint const &_point;
  std::vector<size_t> _dof;
  Eigen::MatrixXd _parametricMesh;
  mutable Eigen::MatrixXd _elementMesh;
  mutable Eigen::MatrixXd _coordinates;
};

inline std::ostream &operator<<(std::ostream &os, PointElementMapper const &mapper)
{
  os << "\nPoint Element Mapper:" << std::endl;
  os << "*** parametric mesh:\n" << mapper.parametricMesh() << std::endl;
  os << "*** element mesh:\n" << mapper.elementMesh() << std::endl;
  return os;
}
