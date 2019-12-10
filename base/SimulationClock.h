#pragma once

#include "Singleton.h"

class SimulationClock : public Singleton<SimulationClock>
{
public:
  void set(double t)
  {
    _last = _time; 
    _time = t;
    _step = 0;
  }

  double time() const {
    return _time;
  }

  double increment() const { return _time - _last; }

  void operator+=(double t) {
    _time += t;
    _step++;
  }

  size_t step() const { return _step; }

private:
  double _last;
  double _time;
  size_t _step;
};
