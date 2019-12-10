#pragma once

class VariableUpdater
{
public:
  virtual ~VariableUpdater() = default;
  virtual bool update() const = 0;
  virtual int priority() const { return -1; }
};
