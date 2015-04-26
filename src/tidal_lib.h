#ifndef TIDAL_LIB_H
#define TIDAL_LIB_H
#include "gravitational_consistency.h"

extern int64_t tidal_extra_range;
void free_tidal_trees(void);
void calc_tidal_forces(struct halo_stash *h, double a1, double a2);

#endif /* TIDAL_LIB_H */
