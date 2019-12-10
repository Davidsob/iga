#ifndef NamedObject_h
#define NamedObject_h


#include <string>

class NamedObject
{
public:
    
  virtual ~NamedObject() = default;

  std::string const &getName() const { return _name; }

protected:

  explicit NamedObject(std::string const &name) 
    : _name(name) {}

private:
  std::string const _name;
};

#endif //NamedObject_h