#include "QuadratureMesh.h"

#include "splines/BSpline.h"
#include "splines/Nurbs.h"

#include "IntegrationPoints.h"
#include "ElementMappers.h"
#include "ParametricMesh.h"

namespace
{
  template<typename Shape>
  void getIntegrationPoints(Shape const &shape,
    std::vector<std::vector<IntegrationPoint>> &ipts,
    std::vector<std::unique_ptr<ElementMapperBase>> &mappers);

  template<>
  void getIntegrationPoints
  ( 
    BSplinePoint const &shape,
    std::vector<std::vector<IntegrationPoint>> &ipts,
    std::vector<std::unique_ptr<ElementMapperBase>> &mappers
  )
  {
    using mapper_t = PointElementMapper;

    mappers.push_back(std::make_unique<mapper_t>(shape));
    mappers.back()->updateElementMesh();
    ipts.push_back(iga::integrationPoints(*mappers.back()));
    for (auto &p : ipts.back()) p.update(); // computes jacobian
  }

  template<>
  void getIntegrationPoints
  ( 
    NurbsCurve const &shape,
    std::vector<std::vector<IntegrationPoint>> &ipts,
    std::vector<std::unique_ptr<ElementMapperBase>> &mappers
  )
  {
    using mapper_t = CurveElementMapper;

    auto elu = iga::meshFromSpan(shape.knot).size()-1;
    for (size_t i = 0; i < elu; i++)
    {
      mappers.push_back(std::make_unique<mapper_t>(shape));
      mappers.back()->updateElementMesh(i);
      ipts.push_back(iga::integrationPoints(*mappers.back(),shape.p));
      for (auto &p : ipts.back()) p.update(); // computes jacobian
    }
  }

  template<>
  void getIntegrationPoints
  (
    NurbsSurface const &shape,
    std::vector<std::vector<IntegrationPoint>> &ipts,
    std::vector<std::unique_ptr<ElementMapperBase>> &mappers
  )
  {
    using mapper_t = ManifoldElementMapper;

    auto elu = iga::meshFromSpan(shape.uknot).size()-1;
    auto elv = iga::meshFromSpan(shape.vknot).size()-1;
    for (size_t j = 0; j < elv; j++)
    {
      for (size_t i = 0; i < elu; i++)
      {
        mappers.push_back(std::make_unique<mapper_t>(shape));
        mappers.back()->updateElementMesh(i,j);
        ipts.push_back(iga::integrationPoints(*mappers.back(),shape.p,shape.q));
        for (auto &p : ipts.back()) p.update(); // computes jacobian
      }
    }
  }

  template<>
  void getIntegrationPoints
  (
    NurbsSolid const &shape,
    std::vector<std::vector<IntegrationPoint>> &ipts,
    std::vector<std::unique_ptr<ElementMapperBase>> &mappers
  )
  {
    using mapper_t = SolidElementMapper;

    auto elu = iga::meshFromSpan(shape.uknot).size()-1;
    auto elv = iga::meshFromSpan(shape.vknot).size()-1;
    auto elw = iga::meshFromSpan(shape.wknot).size()-1;

    for (size_t k = 0; k < elw; k++)
    {
      for (size_t j = 0; j < elv; j++)
      {
        for (size_t i = 0; i < elu; i++)
        {
          mappers.push_back(std::make_unique<mapper_t>(shape));
          mappers.back()->updateElementMesh(i,j,k);
          ipts.push_back(iga::integrationPoints(*mappers.back(),shape.p,shape.q,shape.r));
          for (auto &p : ipts.back()) p.update(); // computes jacobian
        }
      }
    }
  }

  void integrationPointsForShape
  (
    GeometricObject const *_obj,
    std::vector<std::vector<IntegrationPoint>> &ipts,
    std::vector<std::unique_ptr<ElementMapperBase>> &mappers
  )
  {
    if (auto shape = dynamic_cast<NurbsCurve const *>(_obj))
    {
      getIntegrationPoints(*shape,ipts,mappers);
    }
    else if (auto shape = dynamic_cast<NurbsSurface const *>(_obj))
    {
      getIntegrationPoints(*shape,ipts,mappers);
    }
    else if (auto shape = dynamic_cast<NurbsSolid const *>(_obj))
    {
      getIntegrationPoints(*shape,ipts,mappers);
    } else if (auto shape = dynamic_cast<BSplinePoint const *>(_obj))
    {
      getIntegrationPoints(*shape,ipts,mappers);
    }
  }
}


void
QuadratureMesh::
initialize() const
{
  if (!_integrationPoints.empty() && !_mappers.empty()) return;
  if (!_integrationPoints.empty()) _integrationPoints.clear();
  if (!_mappers.empty()) _mappers.clear();
  ::integrationPointsForShape(_obj, _integrationPoints, _mappers);
}