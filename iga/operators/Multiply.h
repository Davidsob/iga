#pragma once

struct Multiply
{
  template<typename A, typename B>
  constexpr decltype(auto) operator()(A const &a, B const &b) const
  {
    return a*b;
  }
};
