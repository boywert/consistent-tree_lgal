#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "../src/gravitational_consistency.h"
#include "../src/gravitational_consistency_vars.h"
#include "../src/check_syscalls.h"
#include "../src/halo_io.h"

#define FAST3TREE_TYPE struct tree_halo
#include "fast3tree.c"

double H=0.7;

#define NUM_TYPES 9
#define T_ALL 0
#define T_CEN 1
#define T_SAT 2
#define T_ALL_PHANT 3
#define T_CEN_PHANT 4
#define T_SAT_PHANT 5
#define T_ALL_NOPHANT 6
#define T_CEN_NOPHANT 7
#define T_SAT_NOPHANT 8
#define BPDEX 4

struct halo_stash now = {0};
float max_mvir=0, min_mvir=0, box_size=0;
float max_vmax=0, min_vmax=0;
int64_t custom_box_size = 0;
int64_t num_bins;
int64_t *counts[NUM_TYPES] = {0};
double *sum_logs = NULL;

static inline void add_to_bin(struct tree_halo *h, double logq, int64_t bin) {
  if (bin < 0) return;
  if (bin >= num_bins) return;
  sum_logs[bin]+=logq;
  counts[T_ALL][bin]++;
  if (h->phantom) counts[T_ALL_PHANT][bin]++;
  else counts[T_ALL_NOPHANT][bin]++;
  if (h->pid < 0) { //Central
    counts[T_CEN][bin]++;
    if (h->phantom) counts[T_CEN_PHANT][bin]++;
    else counts[T_CEN_NOPHANT][bin]++;
  }
  else {
    counts[T_SAT][bin]++;
    if (h->phantom) counts[T_SAT_PHANT][bin]++;
    else counts[T_SAT_NOPHANT][bin]++;
  }
}

void print_function(char *filename, int type, float min) {
  FILE *output;
  int64_t i, j;
  float volume_mul = pow(box_size/H, -3)*BPDEX;
  double avg;
  double bpdex = BPDEX;

  if (type == 1) bpdex *= 2.0;
  output = check_fopen(filename, "w");
  if (type == 0) fprintf(output, "# Mass bin avg. (Msun) ");
  else fprintf(output, "# Vmax bin avg. (km/s) ");
  fprintf(output, "all cen sat all_phant cen_phant sat_phant all_nophant cen_nophant sat_nophant");
  if (type == 0) fprintf(output, " Mass bin center (Msun)\n");
  else fprintf(output, " Vmax bin center (km/s)\n");
  for (i=0; i<num_bins; i++) {
    if (counts[T_ALL][i] == 0) avg = pow(10, min + ((float)i+0.5)/(float)bpdex);
    else avg = pow(10, sum_logs[i]/(double)counts[T_ALL][i]);
    fprintf(output, "%.3e ", avg);
    for (j=0; j<NUM_TYPES; j++) {
      if (counts[j][i]==0) fprintf(output, "%.3e ", volume_mul*1e-2);
      else fprintf(output, "%.3e ", counts[j][i]*volume_mul);
    }
    fprintf(output, "%.3e\n", pow(10, min + ((float)i+0.5)/(float)bpdex));
  }
  fclose(output);
}

void calc_mass_function(char *filename) {
  int64_t i, n, skip=0, bin;
  float logmax=16, logmin=10;
  double logm;
  struct tree_halo *h;

  if (max_mvir > 0) logmax = log10(max_mvir/H);
  if (min_mvir > 0) logmin = log10(min_mvir/H);
  logmin = ((float)((int)(logmin*BPDEX)))/((float)BPDEX);
  assert(logmax > logmin);

  num_bins = ((logmax-logmin)*BPDEX + 2);
  n = num_bins*sizeof(int64_t);
  for (i=0; i<NUM_TYPES; i++) {
    counts[i] = check_realloc(counts[i],n, "Allocating mass function bins.");
    memset(counts[i], 0, n);
  }
  sum_logs = check_realloc(sum_logs,sizeof(double)*num_bins,
			   "Allocating mass function bins.");
  memset(sum_logs, 0, sizeof(double)*num_bins);

  for (i=0; i<now.num_halos; i++) {
    h = &(now.halos[i]);
    if (custom_box_size) {
      for (n=0, skip=0; n<3; n++)
	if (h->pos[n] < 0 || h->pos[n] > box_size) skip=1;
      if (skip) continue;
    }

    if (h->mvir <= 0) continue;
    
    logm = log10(h->mvir/H);
    bin = (float)BPDEX*(logm-logmin);
    add_to_bin(h, logm, bin);
  }
  print_function(filename, 0, logmin);
}

void calc_vmax_function(char *filename) {
  int64_t i, n, skip=0, bin;
  float logmax=16, logmin=10;
  double logv;
  struct tree_halo *h;

  for (i=0; i<now.num_halos; i++) {
    h = &(now.halos[i]);
    if (max_vmax < h->vmax) max_vmax = h->vmax;
    if (!min_vmax || min_vmax > h->vmax) min_vmax = h->vmax;
  }
  if (max_vmax > 0) logmax = log10(max_vmax);
  if (min_vmax > 0) logmin = log10(min_vmax);
  logmin = ((float)((int)(logmin*BPDEX*2.0)))/((float)BPDEX*2.0);
  assert(logmax > logmin);

  num_bins = ((logmax-logmin)*BPDEX*2.0 + 2);
  n = num_bins*sizeof(int64_t);
  for (i=0; i<NUM_TYPES; i++) {
    counts[i] = check_realloc(counts[i],n, "Allocating vmax function bins.");
    memset(counts[i], 0, n);
  }
  sum_logs = check_realloc(sum_logs,sizeof(double)*num_bins,
			   "Allocating mass function bins.");
  memset(sum_logs, 0, sizeof(double)*num_bins);

  for (i=0; i<now.num_halos; i++) {
    h = &(now.halos[i]);
    if (custom_box_size) {
      for (n=0, skip=0; n<3; n++)
	if (h->pos[n] < 0 || h->pos[n] > box_size) skip=1;
      if (skip) continue;
    }

    if (h->vmax <= 0) continue;
    
    logv = log10(h->vmax);
    bin = (float)BPDEX*2.0*(logv-logmin);
    add_to_bin(h, logv, bin);
  }

  print_function(filename, 1, logmin);
}

int sort_halo_order(const void *a, const void *b) {
  float c = now.halos[*((int64_t *)a)].vmax;
  float d = now.halos[*((int64_t *)b)].vmax;
  if (c > d) return -1;
  if (d > c) return 1;
  return 0;
}

void find_parents(void) {
  int64_t i, j, *halo_order = NULL;
  struct fast3tree_results *nearest;
  struct fast3tree *halo_tree;
  struct tree_halo *h1, *h2;
  float max_dist = (float)box_size/2.01, range;

  halo_order = check_realloc(halo_order, sizeof(int64_t)*now.num_halos,
			     "Allocating halo order.");
  halo_tree = fast3tree_init(now.num_halos, now.halos);
  for (i=0; i<now.num_halos; i++) {
    now.halos[i].pid = now.halos[i].upid = -1;
    halo_order[i] = i;
  }
  qsort(halo_order, now.num_halos, sizeof(int64_t), sort_halo_order);

  nearest = fast3tree_results_init();
  for (i=0; i<now.num_halos; i++) {
    h1 = &(now.halos[halo_order[i]]);
    if (h1->flags & DEAD_HALO_FLAG) continue;

    IF_PERIODIC {
      range = h1->rvir/1.0e3;
      if (max_dist < range) range = max_dist;
      fast3tree_find_sphere_periodic(halo_tree, nearest, h1->pos, range);
    } else {
      fast3tree_find_sphere(halo_tree, nearest, h1->pos, h1->rvir/1.0e3);
    }

    for (j=0; j<nearest->num_points; j++) {
      h2 = nearest->points[j];
      if (h2->rvir < h1->rvir) {
	h2->pid = h1->id;
	if (h2->upid < 0) h2->upid = h1->id;
      }
    }
  }
  fast3tree_results_free(nearest);
  fast3tree_free(&halo_tree);
  free(halo_order);
}


int main(int argc, char **argv) {
  char buffer[1024];

  if (argc<2) {
    printf("Usage: %s halo_output [box_size] [find_parents] [h]\n", argv[0]);
    exit(1);
  }

  load_halos(argv[1], &now, 1.0, 0);
  if (argc>2) {
    box_size = atoi(argv[2]);
    custom_box_size = 1;
  }

  if (argc>3 && atoi(argv[3])>0) find_parents();
  if (argc>4 && atof(argv[4])>0) H = atof(argv[4]);

  snprintf(buffer, 1024, "%s.mass_function", argv[1]);
  calc_mass_function(buffer);
  snprintf(buffer, 1024, "%s.vmax_function", argv[1]);
  calc_vmax_function(buffer);
  return 0;
}
