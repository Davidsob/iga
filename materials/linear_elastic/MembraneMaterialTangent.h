#pragma once

#include "utils/MatrixTypes.h"

class MembraneMaterialTangent 
{
public:
  using value_t = StaticMatrixR<3,3>;

  explicit MembraneMaterialTangent(double E, double nu)
    : E(E), nu(nu) {}

  ~MembraneMaterialTangent() = default;

  template<typename Point>
  value_t const operator()(Point const &p) const
  {
    size_t const ik[] = {0,1,0};
    size_t const jl[] = {0,1,1};

    ManifoldElementMapper const * mapper =
      dynamic_cast<ManifoldElementMapper const *>(&p.mapper);

    auto const mu  = E/(2.0*(1+nu));
    auto const lam = E*nu/((1+nu)*(1-nu));

    auto const A = mapper->contravariantMetricTensor(p.para);

    value_t D(value_t::Zero());
    for (size_t a = 0; a < 3; a++)
    {
      auto i = ik[a];
      auto j = jl[a];
      for (size_t b = 0; b < 3; b++)
      {
        auto k = ik[b];
        auto l = jl[b];
        D(a,b) = mu*(A(i,k)*A(j,l) + A(i,l)*A(j,k)) + lam*A(i,j)*A(k,l);
      }
    }

    return D;
  }

private:

  double E, nu;
};
