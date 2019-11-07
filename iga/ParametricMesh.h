#pragma once

#include <Eigen/Dense>

#include <vector>

#include "ShapeFunctions.h"
#include "splines/utils/Converters.h"

namespace iga
{

  template<typename Solid>
  Eigen::MatrixXd parametricMesh(Solid const &solid)
  {
    auto const &uknot = solid.uknot;
    auto const &vknot = solid.vknot;
    auto const &wknot = solid.wknot;

    auto const p = solid.p;
    auto const q = solid.q;
    auto const r = solid.r;

    std::vector<std::vector<double>> mesh;
    for (size_t k = r; k < wknot.size()-r; k++)
    {
      for (size_t j = q; j < vknot.size()-q; j++)
      {
        for (size_t i = p; i < uknot.size()-p; i++)
        {
          mesh.push_back({uknot[i], vknot[j], wknot[k]});
        }
      }
    }

    return convert::to<MatrixXd>(mesh);
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

    int const n = solid.uknot.size()-solid.p-1;
    int const m = solid.vknot.size()-solid.q-1;

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