#pragma once

#include <vector>

#include "splines/BSpline.h"
#include "splines/Nurbs.h"

namespace IO
{
  template<typename Matrix>
  inline void spy(Matrix const &A)
  {
    std::cout << "\nSPY ..." <<std::endl;
    for (int i = 0; i < A.rows(); i++)
    {
      for (int j = 0; j < A.cols(); j++)
      {
        std::cout << !algo::equal(A(i,j),0.0); 
      }
      std::cout << std::endl;
    }
  }

  template<typename Shape, typename Solution>
  void geometryWithSolution(Shape const &solid, Solution const &solution, Shape &out)
  {
    out = solid;
    size_t i = 0;
    for (auto &x : out.Q)
      x.push_back(solution[i++]);
  }

  template<typename Solid, typename Solution, typename ...Args>
  inline void writeSolutionToFile(Solid const &solid,
                                  Solution const &solution,
                                  std::string const &file_name,
                                  Args ...args)
  {
    Solid s; geometryWithSolution(solid,solution,s);
    spline_ops::writeToFile(s,file_name,args...);
  }

  inline void writeXYdata(std::vector<double> const &x, std::vector<double> const &y, std::string const &file_name)
  {
    std::ofstream file;
    file.open(file_name);
    file << x.size() << std::endl;
    size_t k = 0;
    for (auto &xi : x) {
      file << xi << " " << y[k++] << std::endl;
    }
    file.close();
  }
}