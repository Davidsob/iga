#pragma once

#include "base/PairedObjectManagerBase.h"
#include "base/StoredVariable.h"
#include "base/Singleton.h"

#include "splines/GeometricObject.h"

#include "BoundaryConditionBase.h"

class BoundaryConditionManager
  : public PairedObjectManagerBase<GeometricObject, BoundaryConditionBase>
  , public Singleton<BoundaryConditionManager>
{
public:
  ~BoundaryConditionManager() = default;

  friend Singleton<BoundaryConditionManager>;

protected:
  explicit BoundaryConditionManager() = default;
};
