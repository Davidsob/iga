#pragma once

#include "WeakForms.h"
#include "NamedObject.h"

class ConstraintBase : public NamedObject
{
public:
  ConstraintBase(std::string const &name) : NamedObject(name) {};
  ~ConstraintBase() = default;

  virtual WeakFormBase const * getResidual() const = 0;
  virtual WeakFormBase const * getJacobian() const = 0;
  
};