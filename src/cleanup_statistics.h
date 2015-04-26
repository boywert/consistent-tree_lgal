#ifndef CLEANUP_STATISTICS_H
#define CLEANUP_STATISTICS_H

#include <stdint.h>
#include "gravitational_consistency.h"

#define NUM_CLEANUP_STATS 13

#define CLEANUP_STATS_TOTAL                    0
#define CLEANUP_STATS_GOOD                     1
#define CLEANUP_STATS_DEAD                     2
#define CLEANUP_STATS_UNLINKED_PHANTOM         3
#define CLEANUP_STATS_MERGER_FLUCTUATION       4
#define CLEANUP_STATS_SHORT_TRACK              5
#define CLEANUP_STATS_TOO_MANY_PHANT           6
#define CLEANUP_STATS_MISTRACKED_SUBHALO       7
#define CLEANUP_STATS_TRACK_REINSTATED         8
#define CLEANUP_STATS_NEEDS_NEW_DESC           9
#define CLEANUP_STATS_NEW_DESC                10
#define CLEANUP_STATS_FORCED_HALO_RETURN      11
#define CLEANUP_STATS_ALREADY_RESUSCITATED    12

void cleanup_clear_stats(void);
void cleanup_print_stats(int64_t out_num);

void count_good_and_dead_halos(struct halo_stash *s);
void log_forced_halo_return(float a1, float a2, struct tree_halo *prog,
			    struct tree_halo *desc);
void log_found_new_descendant(float a1, float a2, struct tree_halo *prog, struct tree_halo *desc);
void log_halo_needs_new_descendant(float a1, struct tree_halo *h);
void log_track_reinstated(float a1, float a2, struct tree_halo *prog, struct tree_halo *desc);
void log_short_track(float a1, struct tree_halo *h);
void log_short_track_special(float a1, struct tree_halo *h);
void log_unlinked_phantom_halo(float a1, struct tree_halo *h);
void log_merger_fluctuation(float a1, struct tree_halo *h);
void log_already_resuscitated(float a1, float a2, struct tree_halo *prog,
			      struct tree_halo *desc);
void log_too_many_phantoms_cleanup(float a1, struct tree_halo *prog);
void log_mistracked_subhalo(float a1, struct tree_halo *h);
#endif /* CLEANUP_STATISTICS_H */
