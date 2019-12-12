#pragma once

#include <algorithm>
#include <list>
#include <typeinfo>
#include <type_traits>
#include <typeindex>
#include <unordered_map>


template<typename Owner,typename Object, class Enable = void>
class PairedObjectManagerBase
{
public:
  using element_t = std::pair<Owner const *, Object *>;
  using iterator = typename std::list<element_t>::iterator;
  using const_iterator = typename std::list<element_t>::const_iterator;

  iterator begin() { return _objects.begin(); } 
  iterator end()   { return _objects.end(); } 

  virtual ~PairedObjectManagerBase()
  {
    for (element_t pair : _objects)
      delete pair.second;
  }
 
  template<typename T>
  T * has(Owner const &owner) const
  {
    auto equal = [&owner](element_t const &o)
    {
      return (o.first == &owner) && dynamic_cast<T*>(o.second);
    };

    auto it1 = _cachedObjects.find(&owner);
    if (it1 != _cachedObjects.end())
    {
      return dynamic_cast<T*>(it1->second);
    }

    auto it2 = std::find_if(_objects.begin(),_objects.end(),equal);

    if (it2 != _objects.end())
    {

      _cachedObjects[&owner] = dynamic_cast<Object*>((*it2).second);
      return dynamic_cast<T*>((*it2).second);
    }

    return nullptr;
  }

  template<typename T, typename ...Args>
  void add(Owner const &owner, Args const&...args) { 
    if (!has<T>(owner))
    {
      _objects.emplace_front(element_t(&owner, dynamic_cast<Object*>(new T(args...))));
    }
  };

  template<typename T>
  bool erase(Owner const &owner)
  {
    auto equal = [&owner](element_t const &o)
    { 
      return (o.first == &owner) && dynamic_cast<T*>(o.second);
    };

    if (has<T>(owner))
    {
      _cachedObjects.erase(&owner); //ptr is the key

      auto it = std::find_if(_objects.begin(),_objects.end(),equal);
      delete (*it).second;
      _objects.erase(it);
      return true;
    }

    return false;
  }

  size_t size() const { return _objects.size(); }

protected:

  size_t cache_size() const { return _cachedObjects.size(); }

  mutable std::list<element_t> _objects;

  using Cache_t = std::unordered_map<Owner const *, Object *>;
  mutable Cache_t _cachedObjects;
};

template<typename Owner, typename Object>
class PairedObjectManagerBase<Owner,Object, typename std::enable_if<std::is_integral<Object>::type>::type >
{};