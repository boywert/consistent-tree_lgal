#ifndef READ_TREE_H
#define READ_TREE_H

#include <stdint.h>
#include "lgal_tree.h"
#ifndef EXTRA_HALO_INFO
#define EXTRA_HALO_INFO
#endif

struct halo {
  float scale;
  int64_t id, num_prog, phantom, pid, upid, mmp;
  int64_t breadth_first_id, depth_first_id, tree_root_id, orig_halo_id, next_coprogenitor_depthfirst_id, last_progenitor_depthfirst_id;
  int32_t snap_num;
  struct halo *desc, *parent, *uparent, *prog, *next_coprog;
  struct halo *nexthalo;
  struct halo *nexthalo_intree;
  int treenr;
  int64_t id_intree;
  int64_t id_infile;
  int64_t descid;
  double accu_mass;
  float  M200c_all;
  float mvir, orig_mvir, rvir, rs, vrms, scale_of_last_MM,
    vmax, pos[3], vel[3], J[3], spin;
  EXTRA_HALO_INFO
};

struct halo_index_key {
  int64_t id;
  int64_t index;
};

struct halo_list {
  struct halo *halos;
  int64_t num_halos;
  float scale;  
  struct halo_index_key *halo_lookup;
};

struct halo_tree {
  struct halo_list *halo_lists;
  int64_t num_lists;
  int64_t *scale_factor_conv;
  int64_t num_scales;
};

struct lgal_halo_tree {
  int64_t num_halos, num_trees;
  int64_t *num_halos_tree;
  struct halo **root;
  struct halo **lastleaf;
};

extern struct halo_tree halo_tree;
extern struct lgal_halo_tree lgal_halo_tree;
extern struct halo_list all_halos;
extern int64_t total_outputs;
extern int64_t *output_numbers;
extern float *output_scales;
struct halo *lookup_halo_in_list(struct halo_list *hl, int64_t id);
struct halo_list *lookup_scale(float scale);
struct halo_list *find_closest_scale(float scale);
void build_lgal_tree();
void output_lgal_tree();
void read_tree(char *filename);
void delete_tree();

#endif /* READ_TREE_H */
