#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <inttypes.h>
#include "universe_time.h"
#include "check_syscalls.h"
#include "cleanup_statistics.h"

extern FILE *logfile;
extern float max_mvir, min_mvir, box_size;
int64_t counts_min_mass_bin=0, counts_max_mass_bin=0, counts_mass_bins=0;
int64_t *counts[NUM_CLEANUP_STATS] = {0};

FILE *metric_output = NULL;

void cleanup_clear_stats(void)
{
  int64_t i;
  counts_min_mass_bin = 
    (min_mvir>0) ? ((int64_t)(log10f(min_mvir))) : 0;
  counts_max_mass_bin = 
    (max_mvir>0) ? ((int64_t)(log10f(max_mvir))) : 0;
  counts_mass_bins = counts_max_mass_bin - counts_min_mass_bin + 1;
  assert(counts_min_mass_bin <= counts_max_mass_bin);

  for (i=0; i<NUM_CLEANUP_STATS; i++) {
    counts[i] = check_realloc(counts[i], sizeof(int64_t)*counts_mass_bins, 
			      "consistency statistics");
    memset(counts[i], 0, sizeof(int64_t)*counts_mass_bins);
  }
}

void cleanup_print_stats(int64_t out_num) {
  FILE *output;
  char buffer[1024];
  int64_t i,j;
  snprintf(buffer, 1024, "%s/cleanup_statistics_%"PRId64".list", OUTBASE, out_num);
  output = check_fopen(buffer, "w");
  fprintf(output, "#T: Total number of halos.\n");
  fprintf(output, "#G: Good halos (with consistent descendants).\n");
  fprintf(output, "#D: Dead halos (without consistent descendants).\n");
  fprintf(output, "#UP: Unlinked phantom halos.\n");
  fprintf(output, "#MF: Merger Fluctuations.\n");
  fprintf(output, "#ST: Short Tracks Truncated (shorter than %"PRId64").\n", MIN_TIMESTEPS_TRACKED);
  fprintf(output, "#TMP: Tracks which were truncated due to the phantom fraction being higher than %f.\n", MAX_PHANTOM_FRACTION);
  fprintf(output, "#MSH: Tracks which were truncated due to a subhalo not being tracked beyond the virial radius of its host.\n");
  fprintf(output, "#TR: Previously severed tracks Reinstated with new MMP.\n");
  fprintf(output, "#NeedsNewD: Halos previously merging into short track.\n");
  fprintf(output, "#NewD: Halos previously merging into short track with new descendants.\n");
  fprintf(output, "#FHR: Forced halo reinstantiation due to isolated halo merger.\n");
  fprintf(output, "#AR: Short tracks already reinstated for isolated halo merger.\n");
  fprintf(output, "#Mass_bin T G D UP MF ST TMP MSH TR NeedsNewD NewD FHR AR\n");

  for (i=0; i<counts_mass_bins; i++) {
    fprintf(output, "%-2"PRId64" ", counts_min_mass_bin+i);
    for (j=0; j<NUM_CLEANUP_STATS; j++)
      fprintf(output, "%7"PRId64" ", counts[j][i]);
    fprintf(output, "\n");
  }
  fclose(output);
}

inline void add_to_mass_bin(int64_t stat, float mass) {
  int64_t bin;
  if (!(mass>0)) return;
  bin = log10f(mass);
  if (bin < counts_min_mass_bin || bin > counts_max_mass_bin) return;
  bin -= counts_min_mass_bin;
  counts[stat][bin]++;
}

void count_good_and_dead_halos(struct halo_stash *s) {
  int64_t i;
  for (i=0; i<s->num_halos; i++) {
    if (s->halos[i].flags & DEAD_HALO_FLAG) add_to_mass_bin(CLEANUP_STATS_DEAD, s->halos[i].mvir);
    else add_to_mass_bin(CLEANUP_STATS_GOOD, s->halos[i].mvir);
    add_to_mass_bin(CLEANUP_STATS_TOTAL, s->halos[i].mvir);
  }
}

void log_too_many_phantoms_cleanup(float a1, struct tree_halo *prog) {
  add_to_mass_bin(CLEANUP_STATS_TOO_MANY_PHANT, prog->mvir);
  fprintf(logfile, "Removed track with too many phantoms (Phantom count: %"PRId64"; tracking count %"PRId64"): (%f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f with descid %"PRId64").\n",
	  prog->num_mmp_phantoms, prog->tracked+1, a1, prog->id, prog->mvir,
	  prog->vmax, prog->np, prog->pos[0], prog->pos[1], prog->pos[2], 
	  prog->vel[0], prog->vel[1], prog->vel[2],
	  prog->descid);
}

void log_forced_halo_return(float a1, float a2, struct tree_halo *prog,
			    struct tree_halo *desc) {
  add_to_mass_bin(CLEANUP_STATS_FORCED_HALO_RETURN, desc->mvir);
  fprintf(logfile, "Halo forced to return: %f %"PRId64" %e -> %f %"PRId64" %e\n", a1, prog->id, prog->mvir, a2, desc->id, desc->mvir);
}

void log_already_resuscitated(float a1, float a2, struct tree_halo *prog,
			    struct tree_halo *desc) {
  add_to_mass_bin(CLEANUP_STATS_ALREADY_RESUSCITATED, prog->mvir);
  fprintf(logfile, "Halo path already resuscitated: %f %"PRId64" %e -> %f %"PRId64" %e\n", a1, prog->id, prog->mvir, a2, desc->id, desc->mvir);
}


void log_found_new_descendant(float a1, float a2, struct tree_halo *prog, struct tree_halo *desc) {
  add_to_mass_bin(CLEANUP_STATS_NEW_DESC, prog->mvir);
  fprintf(logfile, "Found new descendant: %f %"PRId64" %e -> %f %"PRId64" %e\n", a1, prog->id, prog->mvir, a2, desc->id, desc->mvir);
}

void log_halo_needs_new_descendant(float a1, struct tree_halo *h) {
  add_to_mass_bin(CLEANUP_STATS_NEEDS_NEW_DESC, h->mvir);
  fprintf(logfile, "Halo needs new descendant: %f %"PRId64" %e (pid %"PRId64")\n", a1, h->id, h->mvir, h->pid);
}


void log_track_reinstated(float a1, float a2, struct tree_halo *prog, struct tree_halo *desc) {
  add_to_mass_bin(CLEANUP_STATS_TRACK_REINSTATED, desc->mvir);
  fprintf(logfile, "Found new progenitor for previously severed track: %f %"PRId64" %e to %f %"PRId64" %e\n", a1, 
	  prog->id, prog->mvir, a2, desc->id, desc->mvir);
}

void log_short_track(float a1, struct tree_halo *h) {
  add_to_mass_bin(CLEANUP_STATS_SHORT_TRACK, h->mvir);
  fprintf(logfile, "Removed short track: %f %"PRId64" %e (tracked %"PRId64", pid %"PRId64")\n", a1, h->id, h->mvir, h->tracked, h->pid);
}

void log_short_track_special(float a1, struct tree_halo *h) {
  add_to_mass_bin(CLEANUP_STATS_SHORT_TRACK, h->mvir);
  fprintf(logfile, "Removed short track with no clear progenitor: %f %"PRId64" %e (tracked %"PRId64", pid %"PRId64")\n", a1, h->id, h->mvir, h->tracked, h->pid);
}

void log_mistracked_subhalo(float a1, struct tree_halo *h) {
  add_to_mass_bin(CLEANUP_STATS_MISTRACKED_SUBHALO, h->mvir);
  fprintf(logfile, "Removed mistracked sub: %f %"PRId64" %e (tracked %"PRId64", pid %"PRId64")\n", a1, h->id, h->mvir, h->tracked, h->pid);
}

void log_unlinked_phantom_halo(float a1, struct tree_halo *h) {
  add_to_mass_bin(CLEANUP_STATS_UNLINKED_PHANTOM, h->mvir);
  fprintf(logfile, "Removed unlinked phantom halo: %f %"PRId64" %e\n", a1, h->id, h->mvir);
}

void log_merger_fluctuation(float a1, struct tree_halo *h) {
  add_to_mass_bin(CLEANUP_STATS_MERGER_FLUCTUATION, h->mvir);
  fprintf(logfile, "Removed likely merger fluctuation: %f %"PRId64" %e\n", a1, h->id, h->mvir);
}


