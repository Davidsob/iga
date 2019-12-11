#include "GeometryUtils.h"

namespace GeometryUtils
{
  StaticVectorR<2>
  normal2(StaticVectorR<2> const &a, StaticVectorR<2> const &b)
  {
  auto Tau = b - a;
  auto tau = Tau/Tau.norm();
  return StaticVectorR<2>(tau[1], -tau[0]);
  }
}
