#ifndef MpiUtils_h
#define MpiUtils_h

#include <mpi.h>
#include <iterator>

#include "MpiDatatypeTraits.h"
#include "MpiStorageTraits.h"

#include "MatrixTypes.h"

namespace mpm
{

  inline MPI_Comm communicator()
  {
    return MPI_COMM_WORLD;
  }

  inline int getRank()
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
  }

  inline int getProcesses()
  {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
  }

  inline bool isMaster()
  {
    return (getRank() == 0);
  }

  inline bool isSerial()
  {
    return (getProcesses() == 1);
  }

  inline void barrier()
  {
    MPI_Barrier(MPI_COMM_WORLD);
  }

  inline double wallTime()
  {
    return MPI_Wtime();
  }

  template<typename T, class Operator>
  inline void reduce(T const &input, T &output, Operator op, int root)
  {
    MPI_Reduce(&input, &output, 1, mpi_datatype_traits<T>::type(), op, root, communicator());
  }

  template<typename T>
  inline void broadcast(T &input, int root)
  {
    MPI_Bcast(&input, 1, mpi_datatype_traits<T>::type(), root, communicator());
  }

  template<class InputIt, class OutputIt, class Operator>
  inline void reduce(InputIt first, InputIt last, OutputIt out_first, Operator op, int root)
  {
    using T = typename InputIt::value_type;
    auto n = mpi_storage_traits<T>::size(std::distance(first, last));
    auto type = mpi_datatype_traits<T>::type();
    MPI_Reduce(&*first, &*out_first, n, type, op, root, communicator());
  }

  template<class InputIt>
  inline void broadcast(InputIt first, InputIt last, int root)
  {
    using T = typename InputIt::value_type;
    auto n = mpi_storage_traits<T>::size(std::distance(first, last));
    auto type = mpi_datatype_traits<T>::type();
    MPI_Bcast(&*first, n, type, root, communicator());
  }

  template<typename T>
  void custom_sum(void *in, void *inout, int *len, MPI_Datatype * dptr)
  {
    static auto increment = [](void *a) { return static_cast<char *>(a) + sizeof(T); };
    auto n = (*len)/sizeof(T);
    for (int i=0; i < n; i++)
    {
      auto &out = *static_cast<T*>(inout);
      out += *static_cast<T*>(in);
      in = increment(in);
      inout = increment(inout);
    }
  }

  template<typename T, typename Op>
  void reduceAndShare(T const &local, T &global, Op op, int root)
  {
    reduce(local.begin(), local.end(), global.begin(), op, root);
    broadcast(global.begin(), global.end(), root);
  }

  template<typename T>
  struct VectorSum
  {
    explicit VectorSum()
    {
      MPI_Op_create(custom_sum<T>, true, &op);
    }
    MPI_Op op;
  };
}
#endif //MpiUtils_h
