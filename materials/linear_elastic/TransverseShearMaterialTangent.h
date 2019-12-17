#pragma once

#include "utils/MatrixTypes.h"

class TransverseShearMaterialTangent 
{
public:
  using value_t = StaticMatrixR<2,2>;

  TransverseShearMaterialTangent(double E, double nu)
    : E(E), nu(nu) {}

  ~TransverseShearMaterialTangent() = default;

  template<typename Point>
  value_t operator()(Point const &p) const
  {
    double const factor = 5.0/6.0;
    auto const mu  = E/(2.0*(1+nu));
    value_t D; D << mu,0,0,mu;
    
    return factor*D;
  }

private:

  double E, nu;
};
