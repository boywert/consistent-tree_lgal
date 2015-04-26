#ifndef _RPEAK_H_
#define _RPEAK_H_

struct halo_data {
  float p[3];
  float v[3];
  float mass;
  float vmax;
  float mass_peak;
  float vmax_acc;
  float vmax_peak;
  float scale;
  float cumul_nd;
  int64_t id, pid, tree_root_id, depthfirst_id, last_prog_id;
  float rvir, vrms;
  struct halo_data *next;
};

#endif /*_RPEAK_H_*/
