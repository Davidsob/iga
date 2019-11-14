#pragma once

#include <vector>
/**
 * @brief      This class describes a geometric object. Provides abstract interface for basic geometric objects
 */
class GeometricObject 
{
public:
  using Point_t = std::vector<double>;
  virtual ~GeometricObject() = default;

  size_t dim() const { return coordinates().empty() ? 0 : coordinates()[0].size(); }

  virtual std::vector<Point_t> const &coordinates() const = 0;
};


