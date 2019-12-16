#pragma once

class IgaSolverBase
{
public:
  virtual ~IgaSolverBase() = default;
  
  virtual void initializeSolution() = 0;  
  virtual void solve() = 0;
  virtual void finalizeSolution() = 0;  

  void run()
  {
    initializeSolution();
    solve();
    finalizeSolution();
  }
};