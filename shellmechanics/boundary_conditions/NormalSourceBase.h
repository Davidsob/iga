#pragma once

struct NormalSourceBase
{
  virtual StaticVectorR<3> const operator()(double const u, double const v) const = 0;
};
