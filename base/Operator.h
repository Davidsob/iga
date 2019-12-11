#pragma once

template<typename Index,typename T>
class Operator
{
public:
  using value_t = T;
  using index_t = Index;

  class OperatorEngine
  {
  public:
    using value_t = T;
    using index_t = Index;
    virtual ~OperatorEngine() {}
    virtual value_t operator()(index_t const &) const = 0;
  };
};
