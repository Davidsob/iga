#pragma once

#include "NonCopyable.h"

template<class T>
class Singleton
// : public NonCopyable
{
public:

  virtual ~Singleton() = default;

  static T &instance() { 
    static T instance;
    return instance;
  }

protected:
  // Singleton(Singleton const &) = delete;
  // Singleton &operator=(Singleton const &) = delete;
};