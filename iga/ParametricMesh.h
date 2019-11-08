#pragma once

#include <Eigen/Dense>

#include <vector>
#include <algorithm>

#include "splines/Algorithms.h"
#include "splines/utils/Converters.h"

#include <Eigen/Dense>

namespace iga
{

  inline std::vector<double> meshFromSpan(std::vector<double> const &knot)
  {
    std::vector<double> tmp = knot;
    auto it = std::unique(tmp.begin(), tmp.end(), 
      [](double a, double b)
      {
        return algo::equal(a,b);
      }
    );
    tmp.erase(it,tmp.end());
    return tmp; 
  }


  template<typename Solid>
  Eigen::MatrixXd parametricMesh(Solid const &solid)
  {
    std::vector<std::vector<double>> mesh;
    for (auto w : meshFromSpan(solid.wknot))
      for (auto v : meshFromSpan(solid.vknot))
        for (auto u : meshFromSpan(solid.uknot))
          mesh.push_back({u,v,w});

    return convert::to<Eigen::MatrixXd>(mesh);
  }

  template<typename Solid>
  Eigen::MatrixXd parametricElementMesh(size_t i,size_t j,size_t k,
                                        Solid const &solid,
                                        Eigen::MatrixXd const &pmesh)
  {
    std::vector<std::vector<size_t>> const topology
    {
      {i  ,j  ,k},
      {i+1,j  ,k},
      {i+1,j+1,k},
      {i  ,j+1,k},
      {i  ,j  ,k+1},
      {i+1,j  ,k+1},
      {i+1,j+1,k+1},
      {i  ,j+1,k+1},
    };

    auto const n = meshFromSpan(solid.uknot).size();
    auto const m = meshFromSpan(solid.vknot).size();

    auto idx = [n,m](size_t i, size_t j, size_t k)
    {
      return i + j*n + k*(n*m);
    };

    size_t kk = 0;
    Eigen::MatrixXd emesh(8,3);
    for (auto const &p : topology)
    {
      emesh.row(kk++) = pmesh.row(idx(p[0],p[1],p[2]));
    }

    return emesh;
  }
}