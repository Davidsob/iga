#pragma once

#include "BaseObjectManager.h"
#include "VariableUpdater.h"
#include "Singleton.h"

class VariableUpdateManager
  : public BaseObjectManager<VariableUpdater>
  , public Singleton<VariableUpdateManager>
{
public:
  ~VariableUpdateManager() = default;

  bool update()
  {
    bool success = true;
    if (_sorted.empty()) get_sorted(_sorted);
    for (auto &updater: _sorted)
    {
      success &= updater->update();
    }
    return success;
  }

  friend Singleton<VariableUpdateManager>;

protected:

  void get_sorted(std::list<VariableUpdater *> &sorted)
  {
    sorted.clear();
    sorted = this->_objects;
    sorted.sort([](VariableUpdater const *a, VariableUpdater const *b)
      {
        return a->priority() < b->priority();
      }
    );
  }
  std::list<VariableUpdater*> _sorted;
  explicit VariableUpdateManager() = default;
};