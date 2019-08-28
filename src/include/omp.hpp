#ifndef INCLUDE_OMP_HPP_
#define INCLUDE_OMP_HPP_

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0; }
inline omp_int_t omp_get_max_threads() { return 1; }
inline omp_int_t omp_get_num_procs() { return 1; }
inline bool omp_in_parallel() { return false; }
inline void omp_set_num_threads(int i) {
   (void)i;
   return;
}
#endif

#endif
