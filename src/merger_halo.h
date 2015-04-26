#ifndef _MERGER_HALO_H_
#define _MERGER_HALO_H_
#include "tree_halo.h"

struct merger_halo {
  int64_t id, descid, mmp, pid, upid, desc_pid, orig_id;
  double scale, desc_scale;
  struct merger_halo *desc, *mmp_halo;
  int64_t phantom;
  int64_t first_in_group, next_in_group, first_coprogenitor, next_coprogenitor;
  float pos[3], vel[3], J[3], spin;
  double mvir, rvir, vmax, vrms, rs, orig_mvir;
  double next_mass, prev_mass, incoming_mass;
  float last_mm;
  float vmax_peak, vmax_acc, mvir_peak, mvir_acc;
  int64_t np, num_prog;
  int64_t breadthfirst_order, depthfirst_order, treeroot_id;
  int64_t last_progenitor_df, next_coprogenitor_df, snapnum;
  double extra_params[MAX_EXTRA_PARAMS];
};

#endif /* _MERGER_HALO_H_ */
