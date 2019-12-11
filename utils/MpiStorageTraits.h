#include <mpi.h>

namespace mpm
{
  template<typename T>
  struct mpi_storage_traits
  { 
    static int size(int n) { return n; }
  };

  template<>
  struct mpi_storage_traits<Eigen::Vector2d>
  { 
    static int size(int n) { return n*sizeof(Eigen::Vector2d); }
  };
}