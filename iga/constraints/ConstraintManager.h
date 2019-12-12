#pragma once

#include "base/PairedObjectManagerBase.h"
#include "base/StoredVariable.h"
#include "base/Singleton.h"

#include "splines/GeometricObject.h"

#include "ConstraintBase.h"

class ConstraintManager
  : public PairedObjectManagerBase<GeometricObject, ConstraintBase>
  , public Singleton<ConstraintManager>
{
public:
  ~ConstraintManager() = default;

  friend Singleton<ConstraintManager>;

protected:
  explicit ConstraintManager() = default;
};
