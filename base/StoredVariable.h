#pragma once

#include "NamedObject.h"

class StoredVariableBase
  : public NamedObject
{

public:

  ~StoredVariableBase() = default; 
  virtual bool initialize() = 0;
  virtual bool isInitialized() const = 0;
  virtual bool clear() = 0;

protected:
  explicit StoredVariableBase(std::string const &name)
    : NamedObject(name)
  {} 
};

template<typename T>
class StoredVariable : public StoredVariableBase
{
public:
  using value_t = T;

  ~StoredVariable() = default;

  bool isInitialized() const override
  {
    return !_data.empty() && _data.size() == _size;
  }

  bool initialize() override
  {
    if (isInitialized()) return true;

    _data.resize(_size, _x0);
    return !_data.empty();
  }

  bool clear() override
  {
    _data.clear();
    return _data.empty();
  }

  size_t size() const { return _data.size(); }
  
  std::vector<T> const &data() const
  {
    return _data;
  }

  std::vector<T> &data()
  {
    return _data;
  }

protected:

  StoredVariable(size_t i, T const x0, std::string const &name)
    : StoredVariableBase(name), _size(i), _x0(x0)
  {
    _data.reserve(_size);
  }

private:
  size_t _size;
  T const _x0;
  std::vector<T> _data;
};

template<typename T>
class StoredLocalPointVariable : public StoredVariableBase
{
public:
  using value_t = T;

  ~StoredLocalPointVariable() = default;

  bool isInitialized() const override
  {
    return !_data.empty() && _data.size() == _size;
  }

  bool initialize() override
  { 
    if (isInitialized()) return true;

    std::vector<value_t> zero(_lp_size, _x0);
    _data.resize(_size,  zero);
    return !_data.empty();
  }

  bool clear() override
  {
    _data.clear();
    return _data.empty();
  }

  size_t size() const { return _data.size(); }
  
  std::vector<std::vector<T>> const &data() const
  {
    return _data;
  }

  std::vector<std::vector<T>> &data()
  {
    return _data;
  }

protected:

  StoredLocalPointVariable(size_t i, size_t lp, T const x0, std::string const &name)
    : StoredVariableBase(name), _size(i), _lp_size(lp), _x0(x0)
  {
    _data.reserve(_size);
  }

private:
  size_t _size;
  size_t _lp_size;
  T const _x0;
  std::vector<std::vector<value_t>> _data;
};