#pragma once

#include "ObjectManagerBase.h"
#include "WeakForms.h"
#include "WeakFormTags.h"
#include "Singleton.h"

class WeakFormManager
  : public ObjectManagerBase<WeakFormBase>
  , public Singleton<WeakFormManager>
{
public:

  ~WeakFormManager() = default;
  friend Singleton<WeakFormManager>;

  void initialize()
  {
    static bool initialized = false;
    if (!initialized)
    {
      auto f    = [](WeakFormBase *form) { return dynamic_cast<LinearFormTag*>(form); };
      auto mass = [](WeakFormBase *form) { return !dynamic_cast<ContributesToInteratia*>(form); };
      _linear_it = std::partition(this->begin(), this->end(), f);
      _mass_it = std::partition(_linear_it, this->end(), mass);
      initialized = true;
    }
  }

  struct  LinearForms;
  LinearForms &linearForms() { return _linear_forms; }

  struct  BilinearForms;
  BilinearForms &bilinearForms() { return _bilinear_forms; }

  struct  StiffnessForms;
  StiffnessForms &stiffnessForms() { return _stiffness_forms; }

  struct  MassForms;
  MassForms &massForms() { return _mass_forms; }

protected:
  explicit WeakFormManager()
    : _linear_forms(*this)
    , _bilinear_forms(*this)
    , _stiffness_forms(*this)
    , _mass_forms(*this)
    {}

  iterator _linear_it;
  iterator _mass_it;

  struct LinearForms
  {
    LinearForms(WeakFormManager &mgr)
      : _mgr(mgr) {}

    virtual ~LinearForms() {};

    iterator begin() { return _mgr.begin(); } 
    iterator end()   { return _mgr._linear_it; }

    WeakFormManager &_mgr;
  } _linear_forms;

  struct BilinearForms
  {
    BilinearForms(WeakFormManager &mgr)
      : _mgr(mgr) {}

    virtual ~BilinearForms() {};

    iterator begin() { return _mgr._linear_it; } 
    iterator end()   { return _mgr.end(); }

    WeakFormManager &_mgr;
  } _bilinear_forms;

  struct StiffnessForms
  {
    StiffnessForms(WeakFormManager &mgr)
      : _mgr(mgr) {}

    virtual ~StiffnessForms() {};

    iterator begin() { return _mgr._linear_it; } 
    iterator end()   { return _mgr._mass_it; }

    WeakFormManager &_mgr;
  } _stiffness_forms;

  struct MassForms
  {
    MassForms(WeakFormManager &mgr)
      : _mgr(mgr) {}

    virtual ~MassForms() {};

    iterator begin() { return _mgr._mass_it; } 
    iterator end()   { return _mgr.end(); }

    WeakFormManager &_mgr;
  } _mass_forms;
};