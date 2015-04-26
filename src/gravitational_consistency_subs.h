#ifndef GRAVITATIONAL_CONSISTENCY_SUBS_H
#define GRAVITATIONAL_CONSISTENCY_SUBS_H

#include "gravitational_consistency.h"

int64_t id_to_index(struct halo_stash h, int64_t id);

void read_outputs(float **output_scales, int64_t **outputs, int64_t *num_outputs);
void gc_load_halos(char *filename, struct halo_stash *h, float scale, int64_t orig_data);
void build_nodesc_halo_tree(struct tree_halo *h, int64_t num_h);

float metric(struct tree_halo *h1, struct tree_halo *h2, float sigma_x, float sigma_v, float vmax_avg, float vmax_tol);
float last_ditch_metric(struct tree_halo *h1, struct tree_halo *h2);

void fork_and_print_halos(char *buffer, struct tree_halo *h, int64_t num_halos, int64_t gzip_halos);
void copy_halos_to_evolved(int64_t output_num, int64_t prune_dead_halos);
void calculate_new_output_order(void);
void tag_major_mergers(void);
void calculate_tracking_stats(float a1, float a2);
void initialize_last_timestep(int64_t last_output_num, float scale);
void clear_halo_stash(struct halo_stash *h);

void mark_unphysical_halos(float a);
void break_spurious_links(float a1, float a2);
void build_mmp_list(int64_t stage, float a1, float a2);
int64_t gen_new_phantom_id(struct halo_stash *h, int64_t last_id);
void create_phantom_progenitors(float a1, float a2);
void mark_dead_halos(float a);
void fix_halos_without_descendants(float a1, float a2);
void fix_halos_without_progenitors(float a1, float a2, int64_t stage);
void sort_process_list_by_vmax_evolved();
void calc_metric_stats(int64_t output_num, float a1, float a2);
void print_tidal_forces(int64_t output_num, float a, struct halo_stash *h, struct halo_stash *evolved);
void regen_output_order(int64_t output_num, float scale);
void break_non_mmp_links(float a1, float a2);


#endif /* GRAVITATIONAL_CONSISTENCY_SUBS_H */
