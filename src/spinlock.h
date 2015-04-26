#ifndef _SPINLOCK_H_
#define _SPINLOCK_H_
#include <inttypes.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/
#include "check_syscalls.h"


typedef int16_t spinlock_t;

#define ALLOC_SPINLOCKS(x,n) { (x) = check_realloc(NULL, sizeof(spinlock_t)*(n), "Allocating spinlocks."); memset((x), 0, sizeof(spinlock_t)*(n)); }

static inline int64_t acquire_spinlock(spinlock_t *x, int64_t n) {
  int64_t j=0;
#ifdef _OPENMP
  int64_t k;
  while (1) {
    while (x[n]) for (k=omp_get_thread_num(); k>=0; k--) j++;
#pragma omp atomic
    x[n]++;
    if (x[n]==1) break;
#pragma omp atomic
    x[n]--;
  }
#endif /* _OPENMP */
  return j;
}

static inline void release_spinlock(spinlock_t *x, int64_t n) {
#ifdef _OPENMP
#pragma omp atomic
  x[n]--;
#endif /* _OPENMP */
}

#endif /*  _SPINLOCK_H_ */

  
