#pragma once

#include "VariableManager.h"

#include "utils/MatrixTypeTraits.h"

template<class StoredVar, class Enable = void>
class StoredVariableOperator
{
public:
  using value_t = DynamicVectorR; 
  
  virtual ~StoredVariableOperator() = default;

  template<typename Index>
  DynamicVectorR const operator()(Index const &i) const
  {
    auto const &dof = i.mapper.dof();
    DynamicVectorR var(dof.size());
    auto svar = VariableManager::instance().has<StoredVar>();
    auto const &data = svar->data();
    for (size_t i = 0; i < dof.size(); i++)
    {
      var[i] = data[dof[i]];
    }
    return var; 
  }
};

template<class StoredVar>
class StoredVariableOperator<StoredVar,typename std::enable_if<is_vector_type<typename StoredVar::value_t>::value()>::type>
{
public:
  using value_t = DynamicMatrixR; 
  
  virtual ~StoredVariableOperator() = default;

  template<typename Index>
  DynamicMatrixR const operator()(Index const &i) const
  {
    auto const &dof = i.mapper.dof();
    DynamicMatrixR var(dof.size(), StoredVar::ndof);
    auto svar = VariableManager::instance().has<StoredVar>();
    auto const &data = svar->data();
    for (size_t i = 0; i < dof.size(); i++)
    {
      var.row(i) = data[dof[i]];
    }
    return var; 
  }
};

template<class StoredLpVar>
class StoredLocalPointVariableOperator
{
public:
  using value_t = typename StoredLpVar::value_t;

  virtual ~StoredLocalPointVariableOperator() = default;

  template<typename Index>
  value_t const&operator()(Index const &i) const
  {
    assert(false);
    return value_t();
    // auto svar = VariableManager::instance().has<StoredLpVar>();
    // auto const &data = svar->data();
    // return data[i.cidx][i.lidx];
  }

  template<typename Index>
  value_t &operator()(Index const &i)
  {
    assert(false);
    return value_t(0);
    // auto svar = VariableManager::instance().has<StoredLpVar>();
    // auto &data = svar->data();
    // return data[i.cidx][i.lidx];
  }
};

// template<class StoredCellVar>
// class StoredCellVariableOperator
// {
// public:
//   using value_t = typename StoredCellVar::value_t;

//   virtual ~StoredCellVariableOperator() = default;

//   template<typename Index>
//   value_t const&operator()(Index const &i) const
//   {
//     auto svar = VariableManager::instance().has<StoredCellVar>();
//     auto const &data = svar->data();
//     return data[i];
//   }

//   template<typename Index>
//   value_t &operator()(Index const &i)
//   {
//     auto svar = VariableManager::instance().has<StoredCellVar>();
//     auto &data = svar->data();
//     return data[i];
//   }
// };

