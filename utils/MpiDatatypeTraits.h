#include <mpi.h>

namespace mpm
{
  template<typename T>
  struct mpi_datatype_traits
  {
    static MPI_Datatype type() { return MPI_BYTE; }
  };

  template<> struct mpi_datatype_traits<int>
  {
    static MPI_Datatype type() { return MPI_INT; }
  };

  template<> struct mpi_datatype_traits<float>
  {
    static MPI_Datatype type() { return MPI_FLOAT; }
  };

  template<> struct mpi_datatype_traits<double>
  {
    static MPI_Datatype type() { return MPI_DOUBLE; }
  };
}