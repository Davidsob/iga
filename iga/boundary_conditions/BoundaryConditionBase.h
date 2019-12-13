#pragma once

#include "base/NamedObject.h"

#include "weakforms/WeakForms.h"

class BoundaryConditionBase : public NamedObject
{
public:
  BoundaryConditionBase(std::string const &name) : NamedObject(name) {};
  ~BoundaryConditionBase() = default;

  virtual WeakFormBase const * getRhs() const = 0;
  virtual WeakFormBase const * getLhs() const { return nullptr; };
  
};