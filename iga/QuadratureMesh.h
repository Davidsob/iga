#pragma once

#include "splines/GeometricObject.h"

#include <vector>
#include <memory>

class IntegrationPoint;
class ElementMapperBase;

class QuadratureMesh
{
public:
  using integration_points = std::vector<IntegrationPoint>;
  using const_iterator = std::vector<integration_points>::const_iterator;

  QuadratureMesh(GeometricObject const *obj)
    : _obj(obj)
  {
  }

  ~QuadratureMesh() = default;
  QuadratureMesh(QuadratureMesh const &other) = delete;
  QuadratureMesh(QuadratureMesh const &&other) = delete;

  const_iterator begin() const {
    initialize();
    return _integrationPoints.begin();
  }
  const_iterator end() const { return _integrationPoints.end(); }

private:
  void initialize() const;

  GeometricObject const *_obj;

  mutable std::vector<integration_points> _integrationPoints;
  mutable std::vector<std::unique_ptr<ElementMapperBase>> _mappers;
};