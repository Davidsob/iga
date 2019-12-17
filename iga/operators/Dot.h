#pragma once

struct Dot
{
  template<typename A, typename B>
  decltype(auto) const operator()(A const &a, B const &b) const
  {
    return a.transpose()*b;
  }
};


