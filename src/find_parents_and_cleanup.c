#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <inttypes.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include "gravitational_consistency.h"
#include "gravitational_consistency_vars.h"
#include "find_parents_and_cleanup.h"
#include "cleanup_statistics.h"
#include "grav_config.h"
#include "distance.h"
#include "universe_time.h"
#include "check_syscalls.h"
#include "halo_io.h"
#include "version.h"
#include "litehash.h"

#define FAST3TREE_PREFIX FPAC
#define FAST3TREE_TYPE struct tree_halo
#include "fast3tree.c"

struct halo_stash now={0}, next={0};
struct fast3tree *halo_tree = NULL;
int64_t *halo_order = NULL;
FILE *logfile = NULL;
struct litehash *lh_ids = NULL, *lh_descids = NULL;
int64_t *ids = NULL, *descids = NULL;
int64_t current_id = 1;

float box_size = 0; /* Automatically set; in comoving Mpc/h */
float max_mvir = 0; /* Automatically set; in Msun */
float min_mvir = 0;

int main(int argc, char **argv) {
  int64_t i, j, num_outputs=0, stage;
  float *output_scales=NULL;
  int64_t *outputs=NULL;
  int stat_loc;
  pid_t pid;
  float a1, a2;
  char buffer[1024];

  if (argc==1) {
    fprintf(stderr, "Consistent Trees, Version %s\n", TREE_VERSION);
    fprintf(stderr, "%s.  See the LICENSE file for redistribution details.\n", TREE_COPYRIGHT);
    fprintf(stderr, "Usage: %s options.cfg\n", argv[0]); exit(1);
  }
  grav_config(argv[1], 0);
  print_timing(_FPAC, NULL);

  init_cosmology(Om, Ol, h0);
  init_time_table(Om, h0);
  halo_tree = fast3tree_init(0, NULL);
  read_outputs(&output_scales, &outputs, &num_outputs);

  clear_halo_stash(&next);
  timed_output(_FPAC, "Loading snapshot %"PRId64"...", outputs[0]);
  snprintf(buffer, 1024, "%s/consistent_%"PRId64".list", OUTBASE, outputs[0]);
  load_halos(buffer, &next, output_scales[0], 0);
  snprintf(buffer, 1024, "%s/dead_%"PRId64".list", OUTBASE, outputs[0]);
  load_halos(buffer, &next, output_scales[0], 1);
  build_id_conv_list(&next);

  for (i=0; i<num_outputs; i++)
  {
    snprintf(buffer, 1024, "%s/cleanup_logfile_%"PRId64".list", OUTBASE, outputs[i]);
    logfile = check_fopen(buffer, "w");

    stage = (i==0) ? 1 :
      (i<num_outputs-1) ? 2 : 3;
    a1 = output_scales[i];
    a2 = (stage < 3) ? output_scales[i+1] : a1;
    timed_output(_FPAC, "** Starting work for snapshot %"PRId64" (a=%f)...", outputs[i], a1);

    clear_halo_stash(&now);
    now = next;
    zero_halo_stash(&next);
    if (stage < 3) {
      timed_output(_FPAC, "Loading snapshot %"PRId64" for comparison...", outputs[i+1]);
      snprintf(buffer, 1024, "%s/consistent_%"PRId64".list", OUTBASE, outputs[i+1]);
      load_halos(buffer, &next, a2, 0);
      snprintf(buffer, 1024, "%s/dead_%"PRId64".list", OUTBASE, outputs[i+1]);
      load_halos(buffer, &next, a2, 1);
      build_id_conv_list(&next);
    }

    cleanup_clear_stats();
    timed_output(_FPAC, "Cleaning up excess phantom halos...");
    cleanup_phantoms(a1, a2, stage);


    halo_order = check_realloc(halo_order, sizeof(int64_t)*now.num_halos, "Halo Order");
    
    timed_output(_FPAC, "Building halo tree...");
    for (j=0; j<now.num_halos; j++) halo_order[j] = j;
    fast3tree_rebuild(halo_tree, now.num_halos, now.halos);
    IF_PERIODIC _fast3tree_set_minmax(halo_tree, 0, box_size);
    build_id_conv_list(&now);
    timed_output(_FPAC, "Deciding halo processing order...");
    qsort(halo_order, now.num_halos, sizeof(int64_t), sort_halo_order);
    timed_output(_FPAC, "Finding parents...");
    find_parents();
    timed_output(_FPAC, "Removing halos with short tracks...");
    cleanup_short_tracks(a1, a2, stage, i, num_outputs);
    timed_output(_FPAC, "Refinding parents...");
    find_parents();
    if (stage < 3) {
      timed_output(_FPAC, "Finding new descendants for halos that need them...");
      cleanup_find_new_descendants(a1, a2);
    }
    //cleanup_no_progenitors(a1, stage);
    timed_output(_FPAC, "Refinding parents...");
    find_parents();

    timed_output(_FPAC, "Flagging major mergers...");
    detect_major_mergers(a1);
    if (stage < 3) tag_children_of_centrals(a1);

    fclose(logfile);
    timed_output(_FPAC, "Translating halo IDs...");
    translate_ids(stage);

    timed_output(_FPAC, "Forking and writing halos for snapshot %"PRId64"...", outputs[i]);
    pid = fork();
    if (pid < 1) {
      snprintf(buffer, 1024, "%s/cleanup_logfile_%"PRId64".list", OUTBASE, outputs[i]);
      logfile = check_fopen(buffer, "a");
      if (!pid) clear_halo_stash(&next);
      count_good_and_dead_halos(&now);
      cleanup_print_stats(outputs[i]);
      cleanup_print_halos(outputs[i], stage);
      fclose(logfile);
      if (!pid) clear_halo_stash(&now);
      gzip_file(buffer);
      if (!pid) return 0;
    }

    while (wait4(-1, &stat_loc, WNOHANG, 0)>0);
    timed_output(_FPAC, "** Done with snapshot %"PRId64".", outputs[i]);
  }
  while (wait4(-1, &stat_loc, 0, 0)>0);
  timed_output(_FPAC, "Successfully finished.");
  close_timing_log();
  return 0;
}

void tag_children_of_centrals(float a) {
  for (int64_t j=0; j<next.num_halos; j++)
    next.halos[j].flags -= (next.halos[j].flags & PROG_IS_CENTRAL_FLAG);
  for (int64_t i=0; i<now.num_halos; i++) {
    if (now.halos[i].flags & DEAD_HALO_FLAG) continue;
    int64_t eindex = id_to_index(next, now.halos[i].descid);
    if (eindex < 0) continue;
    if (!(now.halos[i].flags & MMP_FLAG)) {
      if (now.halos[i].pid < 0 && now.halos[i].mvir > next.halos[eindex].mvir*(MAJOR_MERGER / (1.0+MAJOR_MERGER)))
	next.halos[eindex].last_mm = a;
      continue;
    }
    if (next.halos[eindex].last_mm < now.halos[i].last_mm)
      next.halos[eindex].last_mm = now.halos[i].last_mm;
    if (now.halos[i].pid > -1) continue;
    next.halos[eindex].flags |= PROG_IS_CENTRAL_FLAG;
  }
}

void detect_major_mergers(float a) {
  for (int64_t i=0; i<now.num_halos; i++) {
    if (now.halos[i].flags & DEAD_HALO_FLAG) continue;
    int64_t pindex = id_to_index(now, now.halos[i].pid);
    if (!now.halos[i].last_mm) now.halos[i].last_mm = a;
    if (pindex < 0) continue;
    if (!(now.halos[i].flags & PROG_IS_CENTRAL_FLAG)) continue;
    if (now.halos[i].mvir > now.halos[pindex].mvir*(MAJOR_MERGER / (1.0+MAJOR_MERGER))) 
      now.halos[pindex].last_mm = a;
  }
}


int64_t lookup_new_id(int64_t old_id, struct litehash *lh) {
  if (old_id < 0) return -1;
  int64_t new_id = (int64_t)lh_getval(lh, (void *)(&old_id));
  assert(new_id);
  return(new_id);
}

void translate_ids(int64_t stage) {
  if (stage == 1) {
    lh_ids = new_litehash(8);
    ids = check_realloc(ids, sizeof(int64_t)*now.num_halos, "IDs");
    for (int64_t i=0; i<now.num_halos; i++) {
      now.halos[i].orig_id = now.halos[i].id;
      ids[i] = now.halos[i].orig_id;
      lh_setval(lh_ids, (void *)(ids+i), (void *)current_id);
      current_id++;
    }
  } else {
    free(ids);
    ids = descids;
    free_litehash(lh_ids);
    lh_ids = lh_descids;
  }

  descids = check_realloc(NULL, sizeof(int64_t)*next.num_halos, "DescIDs");
  lh_descids = new_litehash(8);

  for (int64_t j=0; j<next.num_halos; j++) {
    next.halos[j].orig_id = next.halos[j].id;
    descids[j] = next.halos[j].orig_id;
    lh_setval(lh_descids, (void *)(descids+j), (void *)current_id);
    current_id++;
  }
  
  for (int64_t i=0; i<now.num_halos; i++) {
    now.halos[i].id = lookup_new_id(now.halos[i].id, lh_ids);
    now.halos[i].pid = lookup_new_id(now.halos[i].pid, lh_ids);
    now.halos[i].upid = lookup_new_id(now.halos[i].upid, lh_ids);
    if (stage < 3)
      now.halos[i].descid = lookup_new_id(now.halos[i].descid, lh_descids);
    else
      now.halos[i].descid = -1;
  }
}


void cleanup_print_halos(int64_t output_num, int64_t stage)
{
  char buffer[1024];
  int64_t i, num_good_halos=0, num_dead_halos=0;
  FILE *good_halos, *dead_halos, *deleted_tracks, *o;
  struct tree_halo *th;

  snprintf(buffer, 1024, "%s/really_consistent_%"PRId64".list", OUTBASE, output_num);
  good_halos = check_fopen(buffer, "w");
  print_halo(good_halos, NULL); //Print header

  snprintf(buffer, 1024, "%s/really_dead_%"PRId64".list", OUTBASE, output_num);
  dead_halos = check_fopen(buffer, "w");
  print_halo(dead_halos, NULL);

  snprintf(buffer, 1024, "%s/really_deleted_tracks_%"PRId64".list", OUTBASE, output_num);
  deleted_tracks = check_fopen(buffer, "w");
  print_halo(deleted_tracks, NULL);

  for (i=0; i<now.num_halos; i++) {
    th = &(now.halos[i]);
    if (th->flags & DEAD_HALO_FLAG) {
      o = dead_halos;
      num_dead_halos++;
      if (th->flags & (MISTRACKED_SUB_FLAG | TOO_MANY_PHANT_FLAG | SHORT_TRACK_FLAG))
	print_halo(deleted_tracks, th);
    } else {
      o = good_halos;
      num_good_halos++;
    }
    print_halo(o, th);
  }

  fclose(good_halos);
  fclose(dead_halos);
  fclose(deleted_tracks);
  if (logfile) fprintf(logfile, "#Good halos: %"PRId64"\n#Dead halos: %"PRId64"\n", num_good_halos, num_dead_halos);
}


void cleanup_phantoms(float a1, float a2, int stage)
{
  //int64_t i, index, j, eindex, mmp_index;
  float f;

  //Remove phantoms from first timestep
  if (stage==1) {
    {for (int64_t i=0; i<now.num_halos; i++) {
      if (!now.halos[i].phantom) continue;
      now.halos[i].flags |= DEAD_HALO_FLAG;
      log_unlinked_phantom_halo(a1, &(now.halos[i]));
      }}
  }

  if (stage>2) return;

  //Mark dead phantoms at next timestep
  {for (int64_t j=0; j<next.num_halos; j++) {
    next.halos[j].mmp_id = -1;
    }}
  
  {for (int64_t i=0; i<now.num_halos; i++) {
    int64_t eindex = id_to_index(next, now.halos[i].descid);
    if (eindex < 0) continue;
    
    if (!(now.halos[i].flags & DEAD_HALO_FLAG)) {
      int64_t mmp_index = id_to_index(now, next.halos[eindex].mmp_id);
      if ((mmp_index < 0) || (now.halos[i].mvir > now.halos[mmp_index].mvir)) {
	next.halos[eindex].mmp_id = now.halos[i].id;
      }
    }
    }}
  
  //Interpolate phantom properties
  {for (int64_t i=0; i<now.num_halos; i++) {
    int64_t eindex = id_to_index(next, now.halos[i].descid);
    if (eindex < 0) continue;
    
    if (next.halos[eindex].flags & DEAD_HALO_FLAG)
      assert(now.halos[i].flags & DEAD_HALO_FLAG);
    
    if ((now.halos[i].flags & DEAD_HALO_FLAG) &&
	(now.halos[i].phantom) && (next.halos[eindex].mmp_id < 0)) {
      if (next.halos[eindex].phantom) {
	next.halos[eindex].flags |= DEAD_HALO_FLAG;
	log_unlinked_phantom_halo(a2, &(next.halos[eindex]));
      }
    }
    
    if (now.halos[i].id != next.halos[eindex].mmp_id) continue;
    
    if (next.halos[eindex].phantom) {
      f = 1.0/((double)(next.halos[eindex].phantom+1));
#define INTERP(x) next.halos[eindex].x =				\
	f*next.halos[eindex].x + (1.0-f)*now.halos[i].x;
#define SINTERP(x) next.halos[eindex].x =				\
	sqrt(f*pow(next.halos[eindex].x,2) + (1.0-f)*pow(now.halos[i].x,2));
#define CINTERP(x) next.halos[eindex].x =				\
	cbrt(f*pow(next.halos[eindex].x,3) + (1.0-f)*pow(now.halos[i].x,3));
      INTERP(mvir);
      CINTERP(rvir);
      CINTERP(vmax);
      CINTERP(vrms);
      INTERP(rs);
      INTERP(np);
      INTERP(J[0]);
      INTERP(J[1]);
      INTERP(J[2]);
      INTERP(spin);
      for (int64_t j=0; j<EXTRA_PARAMS; j++) {
	if (EXTRA_PARAM_INTERPOLATIONS[j]=='c') { CINTERP(extra_params[j]); }
	else if (EXTRA_PARAM_INTERPOLATIONS[j]=='s') {SINTERP(extra_params[j]);}
	else { INTERP(extra_params[j]); }
      }
#undef INTERP
#undef SINTERP
#undef CINTERP
    }
    }}

  {for (int64_t j=0; j<next.num_halos; j++) {
    if (next.halos[j].flags & MMP_FLAG) continue;
    if (next.halos[j].mmp_id > -1) continue;
    next.halos[j].flags |= DEAD_HALO_FLAG;
    next.halos[j].flags |= MERGER_FLUCTUATION_FLAG;
    log_merger_fluctuation(a2, &(next.halos[j]));
    }}
}

void cleanup_short_tracks(float a1, float a2, int stage, int64_t output_index, int64_t num_outputs)
{
  int64_t no_new_tracks = 0;
  if (PADDING_TIMESTEPS>0) {
    if ((output_index + 1 + MIN_TIMESTEPS_TRACKED) > ((num_outputs - 1) - (PADDING_TIMESTEPS + 1)))
      no_new_tracks = 1;
  } else {
    if ((output_index + 1 + MIN_TIMESTEPS_TRACKED) > (num_outputs - 1)) {
      no_new_tracks = 1;
    }
  }

  if (stage == 1) { //Remove tracks from first stage
    {for (int64_t i=0; i<now.num_halos; i++) {
      if ((now.halos[i].tracked-now.halos[i].num_mmp_phantoms) < MIN_TIMESTEPS_TRACKED) {
	now.halos[i].flags |= (DEAD_HALO_FLAG | SHORT_TRACK_FLAG);
	log_short_track(a1, &(now.halos[i]));
      }
      if ((now.halos[i].tracked+1)*MAX_PHANTOM_FRACTION < now.halos[i].num_mmp_phantoms) {
	now.halos[i].flags |= (DEAD_HALO_FLAG | TOO_MANY_PHANT_FLAG);
	log_too_many_phantoms_cleanup(a1, &(now.halos[i]));
      }
      }}
  }

  if (stage > 2) return;

  //Remove short tracks
  {for (int64_t j=0; j<next.num_halos; j++) {
    if (next.halos[j].flags & DEAD_HALO_FLAG) continue;
    int64_t index = id_to_index(now, next.halos[j].mmp_id);

    if ((next.halos[j].pid > -1) &&
	(((next.halos[j].tracked-next.halos[j].num_mmp_phantoms) < MIN_TIMESTEPS_SUB_TRACKED) ||
	 (next.halos[j].tracked_single_mmp < MIN_TIMESTEPS_SUB_MMP_TRACKED)) &&
	((index < 0) || (now.halos[index].flags & MISTRACKED_SUB_FLAG))) {
      if (!(next.halos[j].flags & DEAD_HALO_FLAG))
	log_mistracked_subhalo(a2, &(next.halos[j]));
      next.halos[j].flags |= (DEAD_HALO_FLAG | MISTRACKED_SUB_FLAG);
      next.halos[j].mmp_id = -1;
    }

    if (((index < 0) || (now.halos[index].mvir*30 < next.halos[j].mvir)) &&
	(next.halos[j].mvir > MASS_RES_OK)) { // && (next.halos[j].tracked_single_mmp < MIN_TIMESTEPS_TRACKED)) { 
      // No sensible progenitor found!!!
      next.halos[j].flags |= SHORT_TRACK_FLAG;
      next.halos[j].flags |= DEAD_HALO_FLAG;
      next.halos[j].mmp_id = -1;
      log_short_track_special(a2, &(next.halos[j]));
    }

    if (((next.halos[j].tracked-next.halos[j].num_mmp_phantoms) < MIN_TIMESTEPS_TRACKED || no_new_tracks) &&
	((index < 0 || really_beyond_mmp_ratio(now.halos[index], next.halos[j])) || 
	 (now.halos[index].flags & (SHORT_TRACK_FLAG | DEAD_HALO_FLAG)))) {
      next.halos[j].flags |= SHORT_TRACK_FLAG;
      next.halos[j].flags |= DEAD_HALO_FLAG;
      next.halos[j].mmp_id = -1;
      log_short_track(a2, &(next.halos[j]));
    }
    if (((next.halos[j].tracked+1)*MAX_PHANTOM_FRACTION < next.halos[j].num_mmp_phantoms) &&
	((index < 0 || really_beyond_mmp_ratio(now.halos[index], next.halos[j])) || 
	 (now.halos[index].flags & (TOO_MANY_PHANT_FLAG | DEAD_HALO_FLAG)))) {
      next.halos[j].flags |= (DEAD_HALO_FLAG | TOO_MANY_PHANT_FLAG);
      log_too_many_phantoms_cleanup(a2, &(next.halos[j]));
    }
    }}

  //Check to see if there's an alternate progenitor which could serve as the MMP
  {for (int64_t i=0; i<now.num_halos; i++) {
    if (now.halos[i].flags & DEAD_HALO_FLAG) continue;
    int64_t eindex = id_to_index(next, now.halos[i].descid);
    if (eindex < 0) continue;
    if (!(next.halos[eindex].flags & (SHORT_TRACK_FLAG | TOO_MANY_PHANT_FLAG | MISTRACKED_SUB_FLAG))) continue;
    if (!really_beyond_mmp_ratio(now.halos[i], next.halos[eindex])) {
      int64_t pindex = id_to_index(now, next.halos[eindex].mmp_id);
      if ((pindex < 0) || (now.halos[i].mvir > now.halos[pindex].mvir)) {
	next.halos[eindex].flags -= (next.halos[eindex].flags & DEAD_HALO_FLAG);
	next.halos[eindex].mmp_id = now.halos[i].id;
      }
    }
    }}

  //Reinstate halos which have alternate progenitors
  {for (int64_t j=0; j<next.num_halos; j++) {
    if (!(next.halos[j].flags & (SHORT_TRACK_FLAG | TOO_MANY_PHANT_FLAG | MISTRACKED_SUB_FLAG))) continue;
    if ((next.halos[j].flags & DEAD_HALO_FLAG)) continue;
    int64_t pindex = id_to_index(now, next.halos[j].mmp_id);
    assert(pindex >= 0);
    next.halos[j].flags -= (next.halos[j].flags & (SHORT_TRACK_FLAG | TOO_MANY_PHANT_FLAG | MISTRACKED_SUB_FLAG));
    log_track_reinstated(a1, a2, &(now.halos[pindex]), &(next.halos[j]));
    }}

  //Break links where descendant halo had a short track
  {for (int64_t i=0; i<now.num_halos; i++) {
    if (now.halos[i].flags & DEAD_HALO_FLAG) continue;
    int64_t eindex = id_to_index(next, now.halos[i].descid);
    if (eindex < 0) continue;
    if (next.halos[eindex].flags & (SHORT_TRACK_FLAG | TOO_MANY_PHANT_FLAG | MISTRACKED_SUB_FLAG)) {
      now.halos[i].flags |= (DEAD_HALO_FLAG | FIND_NEW_DESC_FLAG);
      log_halo_needs_new_descendant(a1, &(now.halos[i]));
    }
    }}
}


void cleanup_find_new_descendants(float a1, float a2)
{
  {for (int64_t j=0; j<next.num_halos; j++)
      next.halos[j].tidal_id = -1;
  }
  
  //Try to link up halos to parents
  {for (int64_t i=0; i<now.num_halos; i++) {
    if (!(now.halos[i].flags & FIND_NEW_DESC_FLAG)) continue;
    int64_t pindex = id_to_index(now, now.halos[i].pid);
    if (pindex < 0) continue;
    int64_t eindex = id_to_index(next, now.halos[pindex].descid);
    if (eindex < 0) continue;
    assert(!(next.halos[eindex].flags & DEAD_HALO_FLAG));

    int64_t oeindex = id_to_index(next, now.halos[i].descid);
    if (oeindex >= 0) {
      int64_t tindex = id_to_index(next, next.halos[oeindex].tidal_id);
      if ((tindex < 0) || (next.halos[tindex].mvir < next.halos[eindex].mvir))
	next.halos[oeindex].tidal_id = next.halos[eindex].id;
    }

    now.halos[i].descid = next.halos[eindex].id;
    now.halos[i].flags -= (now.halos[i].flags & FIND_NEW_DESC_FLAG);
    now.halos[i].flags -= (now.halos[i].flags & DEAD_HALO_FLAG);
    log_found_new_descendant(a1, a2, &(now.halos[i]), &(next.halos[eindex]));
    }}

  //Otherwise, maybe have to reinstate halos :(
  if (next.num_halos) {
    int64_t min_descid, max_descid;
    max_descid = min_descid = next.halos[0].descid;
    for (int64_t j=0; j<next.num_halos; j++) {
      if (next.halos[j].descid < min_descid)
	min_descid = next.halos[j].descid;
      if (next.halos[j].descid > max_descid)
	max_descid = next.halos[j].descid;
    }
    int64_t *massive_by_descid = check_realloc(NULL, sizeof(int64_t)*(max_descid-min_descid+1), "Allocating most-massive descid array.\n");
    for (int64_t j=0; j<(max_descid-min_descid+1); j++) massive_by_descid[j]=-1;
    for (int64_t j=0; j<next.num_halos; j++) {
      if (next.halos[j].flags & DEAD_HALO_FLAG) continue;
      int64_t did = next.halos[j].descid-min_descid;
      int64_t msv = massive_by_descid[did];
      if (msv == -1 || next.halos[j].mvir > next.halos[msv].mvir)
	massive_by_descid[did] = j;
    }

    {for (int64_t i=0; i<now.num_halos; i++) {
	if (!(now.halos[i].flags & FIND_NEW_DESC_FLAG)) continue;
	int64_t eindex = id_to_index(next, now.halos[i].descid);
	assert(eindex>=0);
	now.halos[i].flags -= (now.halos[i].flags & DEAD_HALO_FLAG);
	if ((next.halos[eindex].tidal_id == -1) && (next.halos[eindex].descid > -1)) {
	  int64_t did = next.halos[eindex].descid-min_descid;
	  int64_t msv = massive_by_descid[did];
	  if (msv > -1)
	    next.halos[eindex].tidal_id = next.halos[msv].id;
	}

	int64_t new_eindex = id_to_index(next, next.halos[eindex].tidal_id);
	if (new_eindex >= 0) {
	  now.halos[i].descid = next.halos[eindex].tidal_id;
	  log_found_new_descendant(a1, a2, &(now.halos[i]), &(next.halos[new_eindex]));
	}
	else if (next.halos[eindex].flags & DEAD_HALO_FLAG) {
	  next.halos[eindex].flags -= (next.halos[eindex].flags & DEAD_HALO_FLAG);
	  log_forced_halo_return(a1, a2, &(now.halos[i]), &(next.halos[eindex]));
	} else {
	  log_already_resuscitated(a1, a2, &(now.halos[i]), &(next.halos[eindex]));
	}
      }}
    free(massive_by_descid);
  }

  //Fix MMP IDs
  {for (int64_t i=0; i<now.num_halos; i++) {
    now.halos[i].flags -= (now.halos[i].flags & MMP_FLAG);
    if (now.halos[i].flags & DEAD_HALO_FLAG) continue;
    int64_t eindex = id_to_index(next, now.halos[i].descid);
    if (eindex < 0) continue;
    
    int64_t mmp_index = id_to_index(now, next.halos[eindex].mmp_id);
    if ((mmp_index < 0) || (now.halos[i].mvir > now.halos[mmp_index].mvir) 
	|| (now.halos[mmp_index].flags & DEAD_HALO_FLAG)) {
      next.halos[eindex].mmp_id = now.halos[i].id;
    }
    }}

  {for (int64_t j=0; j<next.num_halos; j++) {
    int64_t index = id_to_index(now, next.halos[j].mmp_id);
    if (index < 0) continue;
    now.halos[index].flags |= MMP_FLAG;
    }}
}


inline int64_t id_to_index(struct halo_stash h, int64_t id) {
  if (id < h.min_id || id > h.max_id) return -1;
  return (h.id_conv[id-h.min_id]);
}

void find_parents(void) {
  int64_t i, j;
  float max_dist = (float)box_size/2.01;

  {
    for (i=0; i<now.num_halos; i++) now.halos[i].pid = now.halos[i].upid = -1;

    struct fast3tree_results *nearest = fast3tree_results_init();
    for (i=0; i<now.num_halos; i++) {
      struct tree_halo *h1 = &(now.halos[halo_order[i]]);
      if (h1->flags & DEAD_HALO_FLAG) continue;
      
      IF_PERIODIC {
	float range = h1->rvir/1.0e3;
	if (max_dist < range) range = max_dist;
	fast3tree_find_sphere_periodic(halo_tree, nearest, h1->pos, range);
      } else {
	fast3tree_find_sphere(halo_tree, nearest, h1->pos, h1->rvir/1.0e3);
      }

      for (j=0; j<nearest->num_points; j++) {
	struct tree_halo *h2 = nearest->points[j];
	if (h2->vmax < h1->vmax) {
	  h2->pid = h1->id;
	  if (h2->upid < 0) h2->upid = h1->id;
	}
      }
    }
    fast3tree_results_free(nearest);
  }
}

void calc_num_prog(int stage) {
  int64_t i, j, eindex;
  if (stage > 2) return;
  for (j=0; j<next.num_halos; j++) next.halos[j].num_prog = 0;
  for (i=0; i<now.num_halos; i++) {
    if (now.halos[i].flags & DEAD_HALO_FLAG) continue;
    eindex = id_to_index(next, now.halos[i].descid);
    if (eindex < 0) continue;
    next.halos[eindex].num_prog++;
  }
}

void cleanup_no_progenitors(float a1, int stage) {
  int64_t i, eindex;
  calc_num_prog(stage);
  if (stage < 2) return;
  for (i=0; i<now.num_halos; i++) {
    if (now.halos[i].flags & DEAD_HALO_FLAG) continue;
    if (now.halos[i].num_prog > 0 || now.halos[i].mvir < MASS_RES_OK) continue;
    eindex = id_to_index(next, now.halos[i].descid);
    if (eindex < 0) continue;
    if (next.halos[eindex].num_prog > 1) {
      now.halos[i].flags |= DEAD_HALO_FLAG;
      log_short_track(a1, &(now.halos[i]));
      next.halos[eindex].num_prog--;
    }
  }
}


int sort_halo_order(const void *a, const void *b) {
  float c = now.halos[*((int64_t *)a)].vmax;
  float d = now.halos[*((int64_t *)b)].vmax;
  if (c > d) return -1;
  if (d > c) return 1;
  return 0;
}

void clear_halo_stash(struct halo_stash *h) {
  if (h->halos) h->halos = (struct tree_halo *)realloc(h->halos, 0);
  if (h->id_conv) h->id_conv = (int64_t *)realloc(h->id_conv, 0);
  h->max_id = h->min_id = h->num_halos = 0;
  //max_mvir = 0;
}

void zero_halo_stash(struct halo_stash *h) {
  h->halos = 0;
  h->id_conv = 0;
  h->max_id = h->min_id = h->num_halos = 0;
  //max_mvir = 0;
}



