#pragma once

#include "BaseObjectManager.h"
#include "StoredVariable.h"
#include "Singleton.h"

class VariableManager
  : public BaseObjectManager<StoredVariableBase>
  , public Singleton<VariableManager>
{
public:
  ~VariableManager() = default;

  bool initialize()
  {
    bool success = true;
    for (auto &var : *this)
    {
      success &= var->initialize();
    }
    return success;
  }

  bool clear()
  {
    bool success = true;
    for (auto &var : *this)
    {
      success &= var->clear();
    }
    return success;
  }

  friend Singleton<VariableManager>;

protected:
  explicit VariableManager() = default;
};