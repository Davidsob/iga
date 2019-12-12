#pragma once


class IterativeSolver
{
public:
  virtual ~IterativeSolver() = default;
  
  virtual size_t currentStep() const { return _step; }; 
  virtual void preIterate() = 0; 
  virtual void iterativeSolve() = 0; 
  virtual void postIterate() = 0; 
  virtual bool converged() const = 0;

  void iterate()
  {
    preIterate();
    iterativeSolve();
    postIterate();
    _step++;
  }

protected:
  IterativeSolver() : _step(0) {}
  size_t _step;
};