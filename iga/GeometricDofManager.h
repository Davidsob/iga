#pragma once

#include "base/Singleton.h"

class GeometricDofManager 
  : public Singleton<GeometricDofManager>
{
public:
  // GeometricDofManager() = default;
  ~GeometricDofManager() = default;

  bool hasShape(GeometricObject const *shape)
  {
    if (!shape ) return false;
    return (std::find(shapes.begin(), shapes.end(),shape) == shapes.end());
  }
  
  bool addShape(GeometricObject const *shape)
  {
    if (!shape ) return false;

    if (std::find(shapes.begin(), shapes.end(),shape) == shapes.end())
    {
      shapes.push_back(shape);
      addNewCtrlPoints();
    } else {
      return false;
    }
    return true;
  }

  std::vector<size_t> const &idsForShape(GeometricObject const *shape) const
  {
    static std::vector<size_t> const zero;
    if (!shape) return zero;

    auto it = std::find(shapes.begin(),shapes.end(),shape);
    if (it != shapes.end())
      return sids[std::distance(shapes.begin(),it)];
    return zero;
  }

  std::vector<size_t> const dofForShape(GeometricObject const *shape, size_t const ndof) const
  {
    std::vector<size_t> dof; dof.reserve(sids.size()*ndof);
    for (auto const id : idsForShape(shape))
    {
      for (size_t i = 0; i < ndof; i++) dof.push_back(ndof*id+i);
    }
    return dof;
  }

  std::vector<GeometricObject const*> shapes; 
  std::vector<size_t> ids;
  std::vector<std::vector<size_t>> sids;
  std::vector<std::vector<double>> ctrlPoints;

private:
  void addNewCtrlPoints()
  {
    std::vector<size_t> tmp;
    for (auto const &x : shapes.back()->coordinates())
    {
      auto it = std::find(ctrlPoints.begin(), ctrlPoints.end(), x);
      if (it == ctrlPoints.end())
      {
        auto id = ctrlPoints.size();
        ids.push_back(id);
        tmp.push_back(id);
        ctrlPoints.push_back(x);
      } else {
        auto id = std::distance(ctrlPoints.begin(),it);
        tmp.push_back(id);
      }
    }

    sids.push_back(tmp);
  }
};
