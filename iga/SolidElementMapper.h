#pragma once

#include "ElementMapperBase.h"

class NurbsSolid;

class SolidElementMapper
 : public ElementMapperBase
{
public:

  explicit SolidElementMapper(NurbsSolid const &solid);

  ~SolidElementMapper() = default;

protected:
  Eigen::RowVectorXd _shape(double x1, double x2, double x3)                  const override;
  Eigen::MatrixXd    _grad(double x1, double x2, double x3)                   const override;
  Eigen::MatrixXd    _Grad(double x1, double x2, double x3)                   const override;
  Eigen::Matrix3d    _parametricJacobian(double x1, double x2, double x3)     const override;
  Eigen::Matrix3d    _physicalJacobian(double x1, double x2, double x3)       const override;
  void               _mapIntegrationPoint(double  x1, double  x2, double  x3,
                                          double &p1, double &p2, double &p3) const override;
  void               _updateElementMesh(size_t i, size_t j, size_t k)               override;
private:
  NurbsSolid const &_solid;
  Eigen::MatrixXd _parametricMesh;
  mutable Eigen::MatrixXd _elementMesh;
};