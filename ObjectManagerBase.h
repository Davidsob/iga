#pragma once

#include <algorithm>
#include <list>
#include <typeinfo>
#include <type_traits>
#include <typeindex>
#include <unordered_map>


template<typename Object, class Enable = void>
class ObjectManagerBase
{
public:

  using iterator = typename std::list<Object *>::iterator;
  using const_iterator = typename std::list<Object *>::const_iterator;

  iterator begin() { return _objects.begin(); } 
  iterator end() { return _objects.end(); } 

  virtual ~ObjectManagerBase()
  {
    for (Object *obj : _objects)
      delete obj;
  }
 
  template<typename T>
  T * has() const
  {
    auto it = _cachedObjects.find(typeid(T));
    if (it != _cachedObjects.end())
    {
      return dynamic_cast<T*>(it->second);
    }

    auto it2 = std::find_if(
                  _objects.begin(),
                  _objects.end(),
                  [](Object *o) { return dynamic_cast<T*>(o); }
                );

    if (it2 != _objects.end())
    {
      _cachedObjects[typeid(T)] = dynamic_cast<Object*>(*it2);
      return dynamic_cast<T*>(*it2);
    }    
    
    return nullptr;
  }

  template<typename T, typename ...Args>
  void add(Args const&...args) { 
    if (!has<T>())
    {
      _objects.emplace_front(dynamic_cast<Object*>(new T(args...)));
    }
  };

  template<typename T>
  bool erase()
  {
    if (has<T>())
    {
      _cachedObjects.erase(typeid(T));

      auto it = std::find_if(
                    _objects.begin(),
                    _objects.end(),
                    [](Object *o) { return dynamic_cast<T*>(o); }
                  );
      delete *it;
      _objects.erase(it);

      return true;
    }
    return false;
  }

  size_t size() const { return _objects.size(); }

protected:

  size_t cache_size() const { return _cachedObjects.size(); }

  using Cache_t = std::unordered_map<std::type_index, Object *>;
  
  mutable std::list<Object *> _objects;
  mutable Cache_t _cachedObjects;
};

template<typename Object>
class ObjectManagerBase<Object, typename std::enable_if<std::is_integral<Object>::type>::type >
{};