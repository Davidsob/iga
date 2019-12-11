#ifndef EnumeratedObject_h
#define EnumeratedObject_h

template<typename T>
class EnumeratedObject
{
public:
  static size_t getId() { return EnumeratedObject<T>::id++; }
private:
  static size_t id;
};

template<typename T> size_t EnumeratedObject<T>::id = 0;

#endif //EnumeratedObject_h