#ifndef TREE_HALO_H
#define TREE_HALO_H

#include <stdint.h>

#define MAX_EXTRA_PARAMS (int64_t)30

struct tree_halo {
  int64_t id, descid, mmp_id, tidal_id, pid, upid, orig_id;
  int64_t flags, phantom;
  float tidal_force, mass_factor;
  float pos[3], vel[3], J[3], spin, last_mm;
  double a[3];
  float mvir, rvir, vmax, vrms, rs;
  float pid_vmax, upid_vmax;
  int64_t np, num_prog;
  int64_t tracked, tracked_single_mmp, num_mmp_phantoms;
  double extra_params[MAX_EXTRA_PARAMS];
};

#endif /*TREE_HALO_H*/
