#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "gravitational_consistency.h"
#include "gravitational_consistency_vars.h"
#include "halo_io.h"
#include "grav_config.h"
#include "check_syscalls.h"
#include "masses.h"
#include "halo_io.h"
#include "resort_outputs.h"
#include "inthash.h"
#ifndef NO_FORK
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>
#endif /* NO_FORK */
#include "version.h"

#define MAX_MEMSIZE 250e9
int64_t MAX_FORKS=6;

struct halo_stash now={0}, evolved = {0}, dead = {0};
float box_size=0;
float max_mvir=0;
float min_mvir=0;
int64_t *sort_order = NULL;
struct inthash *forest=NULL, *forest_evolved=NULL, *forest_now=NULL;
int64_t children = 0;

void wait_for_children(int all) {
  int stat_loc;
  for (; children >= MAX_FORKS; children--)
    if (!wait4(-1, &stat_loc, 0, 0))
      children = 1;
  if (all) {
    while (wait4(-1, &stat_loc, 0, 0)>0);
    children = 0;
  }
}

int main(int argc, char **argv) {
  int64_t i, j, num_outputs=0;
  float *output_scales=NULL;
  int64_t *outputs=NULL;
  char buffer[1024];
  int64_t last_output;

  if (argc==1) {
    fprintf(stderr, "Consistent Trees, Version %s\n", TREE_VERSION);
    fprintf(stderr, "%s.  See the LICENSE file for redistribution details.\n", TREE_COPYRIGHT);
    fprintf(stderr, "Usage: %s options.cfg\n", argv[0]); exit(1);
  }
  if (argc>1) grav_config(argv[1], 1);
  print_timing(_RO, NULL);


  read_outputs(&output_scales, &outputs, &num_outputs);
  gen_ff_cache();
  clear_halo_stash(&now);
  forest = new_inthash();
  forest_evolved = new_inthash();
  forest_now = new_inthash();
  
  for (i = num_outputs-1; i>=0; i--) {
    timed_output(_RO, "** Starting work on snapshot %"PRId64" (a=%f)...", outputs[i], output_scales[i]);
    timed_output(_RO, "Loading snapshot %"PRId64"...", outputs[i]);
    snprintf(buffer, 1024, "%s/really_consistent_%"PRId64".list", OUTBASE, outputs[i]);
    clear_halo_stash(&now);
    load_halos(buffer, &now, output_scales[i], 0);
    MAX_FORKS = (MAX_MEMSIZE / (now.num_halos*sizeof(struct tree_halo)*2.0)) - 1;
    if (MAX_FORKS > 6) MAX_FORKS = 6;
    if (MAX_FORKS < 1) MAX_FORKS = 1;
    if (LIMITED_MEMORY) MAX_FORKS = 1;
    build_id_conv_list(&now);
    if (i==num_outputs - 1) last_output = 1;
    else last_output = 0;
    if (last_output) {
      timed_output(_RO, "Setting up forests...");
      for (j=0; j<now.num_halos; j++) {
	ih_setint64(forest, now.halos[j].id, now.halos[j].id);
	ih_setint64(forest_now, now.halos[j].id, now.halos[j].id);
      }
    } else {
      timed_output(_RO, "Updating forests...");
      free_inthash(forest_evolved);
      forest_evolved = forest_now;
      forest_now = new_inthash();
      for (j=0; j<now.num_halos; j++) {
	int64_t desc_rootid = ih_getint64(forest_evolved, now.halos[j].descid);
	assert(desc_rootid != IH_INVALID);
	ih_setint64(forest_now, now.halos[j].id, desc_rootid);
      }
    }
    
    for (j=0; j<now.num_halos; j++) {
      if (now.halos[j].pid < 0) continue;
      int64_t our_rootid = ih_getint64(forest_now, now.halos[j].id);
      int64_t p_rootid = ih_getint64(forest_now, now.halos[j].pid);
      int64_t up_rootid = ih_getint64(forest_now, now.halos[j].upid);
      assert(up_rootid != IH_INVALID && p_rootid != IH_INVALID 
	     && our_rootid != IH_INVALID);
      link_forests(p_rootid, our_rootid);
      link_forests(up_rootid, our_rootid);
    }

    sort_order = check_realloc(sort_order, sizeof(int64_t)*now.num_halos,
			       "Allocating sort order for halos.");
    for (j=0; j<now.num_halos; j++) sort_order[j] = j;
    timed_output(_RO, "Resorting halos...");
    resort_halos(last_output);
    timed_output(_RO, "Forking and writing halos...");
    snprintf(buffer, 1024, "%s/really_consistent_%"PRId64".list", OUTBASE, outputs[i]);
    fork_and_print_halos(buffer, &now);

    timed_output(_RO, "Loading dead halos...");
    clear_halo_stash(&dead);
    snprintf(buffer, 1024, "%s/really_dead_%"PRId64".list", OUTBASE, outputs[i]);
    load_halos(buffer, &dead, output_scales[i], 1);
    sort_order = check_realloc(sort_order, sizeof(int64_t)*dead.num_halos,
			       "Allocating sort order for halos.");
    for (j=0; j<dead.num_halos; j++) sort_order[j] = j;
    timed_output(_RO, "Resorting dead halos...");
    qsort(sort_order, dead.num_halos, sizeof(int64_t), sort_dead_by_mvir);
    timed_output(_RO, "Forking and writing dead halos...");
    snprintf(buffer, 1024, "%s/really_dead_%"PRId64".list", OUTBASE, outputs[i]);
    fork_and_print_halos(buffer, &dead);
    timed_output(_RO, "** Done with snapshot %"PRId64".", outputs[i]);
  }
  wait_for_children(1);

  timed_output(_RO, "Printing forests...");
  print_forests();
  timed_output(_RO, "Successfully finished.");
  close_timing_log();
  return 0;
}


void print_forests(void) {
  int64_t i;
  int64_t *keys = ih_keylist(forest);
  char buffer[1024];
  snprintf(buffer, 1024, "%s/forests.list", TREE_OUTBASE);
  FILE *output = check_fopen(buffer, "w");
  fprintf(output, "#TreeRootID ForestID\n");
  for (i=0; i<forest->elems; i++)
    fprintf(output, "%"PRId64" %"PRId64"\n", keys[i],
	    find_forest_root(keys[i]));
  fclose(output);
  free(keys);
}

int64_t find_forest_root(int64_t rootid) {
  int64_t newroot = ih_getint64(forest, rootid);
  assert(newroot != IH_INVALID);
  if (newroot != rootid) {
    newroot = find_forest_root(newroot);
    ih_setint64(forest, rootid, newroot);
  }
  return newroot;
}

void link_forests(int64_t rootid1, int64_t rootid2) {
  if (rootid1 == rootid2) return;
  rootid1 = find_forest_root(rootid1);
  rootid2 = find_forest_root(rootid2);
  if (rootid1 == rootid2) return;
  if (rootid1 < rootid2) ih_setint64(forest, rootid2, rootid1);
  else ih_setint64(forest, rootid1, rootid2);
}

void fork_and_print_halos(char *filename, struct halo_stash *h)
{
  char buffer[1024];
  pid_t pid;
  FILE *o;
  int64_t i;
  
  snprintf(buffer, 1024, "%s.new", filename);
  wait_for_children(0);
  pid = fork();
  if (pid < 1) {
    if (!pid) {
      clear_halo_stash(&evolved);
      if (h != &now) clear_halo_stash(&now);
      if (h != &dead) clear_halo_stash(&dead);
    }
    o = check_fopen(buffer, "w");
    print_halo(o, NULL);
    for (i=0; i<h->num_halos; i++) print_halo(o, h->halos + sort_order[i]);      
    fclose(o);
    rename(buffer, filename);
    if (!pid) exit(0);
  } else {
    children++;
  }
}

inline int64_t id_to_index(struct halo_stash h, int64_t id) {
  if (id < h.min_id || id > h.max_id) return -1;
  return (h.id_conv[id-h.min_id]);
}

int sort_by_mvir(const void *a, const void *b) {
  const struct tree_halo *c = now.halos + *((int64_t *)a);
  const struct tree_halo *d = now.halos + *((int64_t *)b);
  if (c->mvir > d->mvir) return -1;
  if (d->mvir > c->mvir) return 1;
  if (c->id < d->id) return -1;
  return 1;
}

int sort_dead_by_mvir(const void *a, const void *b) {
  const struct tree_halo *c = dead.halos + *((int64_t *)a);
  const struct tree_halo *d = dead.halos + *((int64_t *)b);
  if (c->mvir > d->mvir) return -1;
  if (d->mvir > c->mvir) return 1;
  if (c->id < d->id) return -1;
  return 1;
}

int sort_by_desc(const void *a, const void *b) {
  const struct tree_halo *c = now.halos + *((int64_t *)a);
  const struct tree_halo *d = now.halos + *((int64_t *)b);
  int64_t cindex = id_to_index(evolved, c->descid);
  int64_t dindex = id_to_index(evolved, d->descid);
  if (cindex < 0 && dindex >= 0) return 1;
  if (cindex >=0 && dindex < 0) return -1;
  if (cindex < dindex) return -1;
  if (dindex < cindex) return 1;
  if (c->mvir < d->mvir) return 1;
  if (d->mvir < c->mvir) return -1;
  if (c->id < d->id) return -1;
  return 1;
}


void clear_halo_stash(struct halo_stash *h) {
  if (h->halos) h->halos = (struct tree_halo *)realloc(h->halos, 0);
  if (h->id_conv) h->id_conv = (int64_t *)realloc(h->id_conv, 0);
  h->max_id = h->min_id = h->num_halos = 0;
}

void zero_halo_stash(struct halo_stash *h) {
  h->halos = 0;
  h->id_conv = 0;
  h->max_id = h->min_id = h->num_halos = 0;
}

void resort_halos(int last_output)
{
  int64_t i;
  if (last_output)
    qsort(sort_order, now.num_halos, sizeof(int64_t), sort_by_mvir);
  else
    qsort(sort_order, now.num_halos, sizeof(int64_t), sort_by_desc);
  clear_halo_stash(&evolved);
  evolved.halos = check_realloc(evolved.halos, sizeof(struct tree_halo)*now.num_halos,
				"Allocating space for sorted halos.");
  evolved.num_halos = now.num_halos;
  for (i=0; i<now.num_halos; i++) evolved.halos[i] = now.halos[sort_order[i]];
  build_id_conv_list(&evolved);
  return;
}
