#pragma once

#include <vector>

#include "splines/BSpline.h"
#include "splines/Nurbs.h"

namespace IO
{

  template<typename Shape, typename Solution>
  void geometryWithSolution(Shape const &solid, Solution const &solution, Shape &out)
  {
    out = solid;
    size_t i = 0;
    for (auto &x : out.Q)
      x.push_back(solution[i++]);
  }

  template<typename Solid>
  inline void writeSolutionToFile(Solid const &solid,
                                  std::vector<double> const &solution,
                                  std::string const &file_name,
                                  int ulevel = 10, int vlevel=10, int wlevel=10)
  {
    Solid s; geometryWithSolution(solid,solution,s);
    spline_ops::writeToFile(s,file_name,ulevel,vlevel,wlevel);
  }
}