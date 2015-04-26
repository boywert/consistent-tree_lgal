#ifndef _GRAVITATIONAL_STATISTICS_H
#define _GRAVITATIONAL_STATISTICS_H

#include "gravitational_consistency.h"

#define NUM_STATS 17
#define STATS_TOTAL                      0
#define STATS_GOOD                       1
#define STATS_DEAD                       2
#define STATS_NO_DESC                    3
#define STATS_DESC_NOT_FOUND             4
#define STATS_PROG_GRAV_INCONSISTENCY    5
#define STATS_DESC_GRAV_INCONSISTENCY    6
#define STATS_INSUFFICIENT_TIDAL         7
#define STATS_PROG_GRAV_REPAIR           8
#define STATS_DESC_GRAV_REPAIR           9
#define STATS_TIDAL_REPAIR               10
#define STATS_PHANTOM                    11
#define STATS_NEW_PHANTOM                12
#define STATS_REM_NODESC                 13
#define STATS_MASSIVE_TIDAL              14
#define STATS_SPURIOUS_MMP               15
#define STATS_TOO_MANY_PHANT             16


void add_to_mass_bin(int64_t stat, float mass);
void clear_stats(void);
void log_too_many_phantoms(float a1, float a2, struct tree_halo *prog);
void log_phantom_halo(float a1, float a2, struct tree_halo *prog, struct tree_halo *desc);
void log_dead_halo(float a, struct tree_halo *halo);
void log_tidal_repair(float a1, float a2, struct tree_halo *halo, struct tree_halo *target);
void log_grav_repair(float a1, float a2, float min_metric, struct tree_halo *prog, struct tree_halo *desc);
void log_desc_not_found(float a1, float a2, struct tree_halo *prog);
void log_no_desc(float a, struct tree_halo *halo);
void log_not_enough_tidal(float a1, float a2, struct tree_halo *prog, struct tree_halo *desc);
void log_break_massive_tidal(float a1, float a2, struct tree_halo *prog, struct tree_halo *desc);
void log_grav_inconsistency(float a1, float a2, float metric_val, struct tree_halo *prog, struct tree_halo *desc);
void log_spurious_mmp_link_broken(float a1, float a2, struct tree_halo *prog, struct tree_halo *desc);
void print_stats(int64_t out_num);
void count_good_halos(struct halo_stash *s);
void turn_on_full_metric_output(int64_t output_num, float a1, float a2);
void build_metric_stats(struct tree_halo *prog, struct tree_halo *parent, struct tree_halo *desc);
void finish_metric_stats(int64_t output_num, float a1, float a2);
void sigma_x_v_vmax(struct tree_halo *h, float *s_x, float *s_v,
		    float *vmax_avg, float *s_vmax, float a);



#endif /*_GRAVITATIONAL_STATISTICS_H */
