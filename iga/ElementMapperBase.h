#pragma once

#include <Eigen/Dense>

#include <vector>

class GeometricObject; 

class ElementMapperBase
{
public:
  virtual ~ElementMapperBase() = default;

  virtual std::vector<size_t> const &dof() const = 0;
  virtual GeometricObject const * geometry() const = 0;
  virtual Eigen::MatrixXd const &coordinates() const = 0;

  template<typename ParametricPoint>
  Eigen::RowVectorXd shape(ParametricPoint const &p) const
  {
    return _shape(p[0],p[1],p[2]);
  }

  template<typename ParametricPoint>
  Eigen::MatrixXd Grad(ParametricPoint const &p) const
  {
    return _Grad(p[0],p[1],p[2]);
  }

  template<typename ParametricPoint>
  Eigen::MatrixXd grad(ParametricPoint const &p) const
  {
    return _grad(p[0],p[1],p[2]);
  }

  template<typename IsoparametricPoint>
  Eigen::MatrixXd parametricJacobian(IsoparametricPoint const &p) const
  {
    return _parametricJacobian(p[0],p[1],p[2]);
  }

  template<typename ParametricPoint>
  Eigen::MatrixXd physicalJacobian(ParametricPoint const &p) const
  {
    return _physicalJacobian(p[0],p[1],p[2]);
  }

  void updateElementMesh(size_t i=0,size_t j=0,size_t k=0)
  {
    _updateElementMesh(i,j,k);
  }

  template<typename IntegrationPoint>
  void mapIntegrationPoint(IntegrationPoint &p) const
  {
    _mapIntegrationPoint(p.gauss[0],p.gauss[1],p.gauss[2],
                         p.para[0], p.para[1], p.para[2]);

    // set weight of integration point
    auto const jac_iso  = parametricJacobian(p.gauss).determinant();
    auto const jac_phys = physicalJacobian(p.para).determinant();

    // std::printf("jac_iso = %f, jac_phys = %f\n",jac_iso, jac_phys);
    p.jdet = jac_phys*jac_iso;
  }

protected:

  explicit ElementMapperBase() = default;

  virtual Eigen::RowVectorXd _shape(double x1, double x2, double x3) const = 0;
  virtual Eigen::MatrixXd    _grad(double x1, double x2, double x3) const = 0;
  virtual Eigen::MatrixXd    _Grad(double x1, double x2, double x3) const = 0;
  virtual Eigen::MatrixXd    _parametricJacobian(double x1, double x2, double x3) const = 0;
  virtual Eigen::MatrixXd    _physicalJacobian(double x1, double x2, double x3) const = 0;
  virtual void               _mapIntegrationPoint(double  x1, double  x2, double  x3,
                                                  double &p1, double &p2, double &p3) const = 0;
  virtual void               _updateElementMesh(size_t i, size_t j, size_t k) = 0;
};