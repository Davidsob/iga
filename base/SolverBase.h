#pragma once

class SolverBase
{
public:
  virtual ~SolverBase() = default;
  
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