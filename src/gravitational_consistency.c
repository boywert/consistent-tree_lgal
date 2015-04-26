#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <inttypes.h>
#ifndef NO_FORK
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>
#endif /* NO_FORK */
#include "gravitational_consistency.h"
#include "gravitational_consistency_vars.h"
#include "gravitational_consistency_subs.h"
#include "gravitational_statistics.h"
#include "distance.h"
#include "universe_time.h"
#include "grav_config.h"
#include "masses.h"
#include "tidal_lib.h"
#include "halo_evolve_lib.h"
#include "check_syscalls.h"
#include "halo_io.h"
#include "version.h"

#define FAST3TREE_TYPE struct halo_pos
#define FAST3TREE_PREFIX GC
#include "fast3tree.c"

struct halo_stash now={0}, evolved={0};
int64_t *output_order=NULL, *output_order_inv=NULL;
int64_t *process_order=NULL;
struct fast3tree *halo_tree = NULL;
struct halo_pos *halo_positions = NULL;
float box_size=0;
float max_mvir=0;
float min_mvir=0;
FILE *logfile = NULL;

int main(int argc, char **argv) {
  int64_t i, j, num_outputs=0;
  float *output_scales=NULL;
  int64_t *outputs=NULL;
  char buffer[1024];
  float a1, a2;
  int64_t restart_number = 0;
  int stat_loc;
  int64_t full_outputs = 0;


  if (argc==1) {
    fprintf(stderr, "Consistent Trees, Version %s\n", TREE_VERSION);
    fprintf(stderr, "%s.  See the LICENSE file for redistribution details.\n", TREE_COPYRIGHT);
    fprintf(stderr, "Usage: %s options.cfg\n", argv[0]); exit(1);
  }
  if (argc>1) grav_config(argv[1], 1);
  if (argc>2) restart_number = atoi(argv[2]);
  print_timing(_GC, NULL);

  init_cosmology(Om, 1.0-Om, h0);
  init_time_table(Om, h0);
  halo_tree = fast3tree_init(0, NULL);
  read_outputs(&output_scales, &outputs, &num_outputs);
  gen_ff_cache();

  i=num_outputs-1;
  if (restart_number>0) while (i>0 && outputs[i]!=restart_number) i--;

  //Reprint last output step in new format, get halo output order.
  if (i==num_outputs-1) {
    initialize_last_timestep(outputs[i], output_scales[i]);
    i--;
  }
  else {
    regen_output_order(outputs[i+1], output_scales[i+1]);
  }

  for (; i>=0; i--)
  {
    a1 = output_scales[i];
    a2 = output_scales[i+1];
    timed_output(_GC, "** Starting work on snapshot %"PRId64" (a=%f)...", outputs[i], a1);

    //Decide whether to output full info:
    if (!(((num_outputs-2)-i)%10)) tidal_extra_range = full_outputs = 1;
    else full_outputs = tidal_extra_range = 0;

    //Gravitationally evolve halos at the next timestep to their
    // positions at the current timestep
    timed_output(_GC, "Evolving halos to previous timestep...");
    evolve_halos(a2, a1, &evolved);
    build_id_conv_list(&evolved);
    snprintf(buffer, 1024, "%s/evolved_%"PRId64".list", OUTBASE, outputs[i]);    
    fork_and_print_halos(buffer, evolved.halos, evolved.num_halos, !full_outputs);

    //Start up logfile
    snprintf(buffer, 1024, "%s/logfile_%"PRId64".list", OUTBASE, outputs[i]);
    logfile = check_fopen(buffer, "w");
    fprintf(logfile, "#Links from %f to %f\n#Dt: %f Myr\n", a1, a2,
	    fabs(scale_to_years(a1)-scale_to_years(a2)) / 1.0e6);


    //Load halos at current timestep
    max_mvir = min_mvir = 0;
    timed_output(_GC, "Loading snapshot %"PRId64" for comparison...", outputs[i]);
    snprintf(buffer, 1024, "%s/out_%"PRId64".list%s", INBASE, outputs[i],
	     ((!strcasecmp(INPUT_FORMAT, "BINARY")) ? ".bin" : ""));
    clear_halo_stash(&now);
    gc_load_halos(buffer, &now, a1, 1);
    timed_output(_GC, "Calculating tidal forces...");
    calc_tidal_forces(&now, a1, a2);
    for (j=0; j<evolved.num_halos; j++) {
      if (evolved.halos[j].mvir > max_mvir) max_mvir = evolved.halos[j].mvir;
      if (evolved.halos[j].mvir < min_mvir) min_mvir = evolved.halos[j].mvir;
    }
    clear_stats();

    //Mark halos which have no mass or vmax as dead.
    timed_output(_GC, "Removing halos with M=0 or Vmax=0...");
    mark_unphysical_halos(a1);

    //Build most-massive progenitor links
    timed_output(_GC, "Building MMP links...");
    build_mmp_list(1, a1, a2);

    //Build metric statistics
    if (full_outputs) {
      timed_output(_GC, "Printing full tidal forces...");
      turn_on_full_metric_output(outputs[i], a1, a2);
      print_tidal_forces(outputs[i], a1, &now, &evolved);
    }

    timed_output(_GC, "Breaking non-MMP links...");
    break_non_mmp_links(a1, a2);
    timed_output(_GC, "Comparing predicted evolution to actual halos...");
    calc_metric_stats(outputs[i], a1, a2);

    //Break links which are gravitationally inconsistent
    timed_output(_GC, "Breaking gravitationally inconsistent links...");
    break_spurious_links(a1, a2);

    //Find links which are gravitationally consistent, if possible
    timed_output(_GC, "Calculating halo processing order...");
    sort_process_list_by_vmax_evolved();
    timed_output(_GC, "Fixing halos without progenitors (stage 1)...");
    fix_halos_without_progenitors(a1, a2, 1);
    timed_output(_GC, "Fixing halos without progenitors (stage 2)...");
    fix_halos_without_progenitors(a1, a2, 2);
    timed_output(_GC, "Fixing halos without progenitors (stage 3)...");
    fix_halos_without_progenitors(a1, a2, 3);
    timed_output(_GC, "Fixing halos without descendants...");
    fix_halos_without_descendants(a1, a2);

    //Mark halos which still have no descendants as dead halos
    timed_output(_GC, "Marking halos with no descendants as dead...");
    mark_dead_halos(a1);

    //Create phantom halos for halos which have no progenitors
    timed_output(_GC, "Creating phantom halos...");
    create_phantom_progenitors(a1, a2);

    //Regenerate most-massive progenitor links, and
    //recheck to make sure that all non-dead halos have descendants.
    timed_output(_GC, "Rebuilding MMP list and tagging major mergers...");
    build_mmp_list(2, a1, a2);
    tag_major_mergers();

    //Remove tracks with too many phantom halos
    timed_output(_GC, "Removing halos with too many phantoms...");
    calculate_tracking_stats(a1, a2);

    timed_output(_GC, "Recalculating output order...");
    calculate_new_output_order();
    timed_output(_GC, "Calculating stats for halos added/removed...");
    count_good_halos(&now);
    print_stats(outputs[i]);
    timed_output(_GC, "Writing revised halo lists...");
    copy_halos_to_evolved(outputs[i], (i < (num_outputs - 
				  PADDING_TIMESTEPS - 1)) ? 1 : 0);
    while (wait4(-1, &stat_loc, WNOHANG, 0)>0);
    timed_output(_GC, "** Done with snapshot %"PRId64".", outputs[i]);
  }
  while (wait4(-1, &stat_loc, 0, 0)>0);
  timed_output(_GC, "Successfully finished.");
  close_timing_log();
  return 0;
}


inline int64_t id_to_index(struct halo_stash h, int64_t id) {
  if (id < h.min_id || id > h.max_id) return -1;
  return (h.id_conv[id-h.min_id]);
}

//Generates a phantom id by filling in holes in the id list if possible
int64_t gen_new_phantom_id(struct halo_stash *h, int64_t last_id) {
  int64_t new_id = last_id + 1;
  if (new_id > h->max_id) return new_id;
  if (new_id < h->min_id) new_id = h->min_id;
  for (; new_id <= h->max_id; new_id++) {
    if (h->id_conv[new_id-h->min_id] < 0) {
      h->id_conv[new_id-h->min_id] = h->num_halos;
      return new_id;
    }
  }
  return(h->max_id + 1);
}

void tag_major_mergers(void) {
  int64_t i;
#pragma omp parallel for private(i) schedule(static,1000)
  for (i=0; i<now.num_halos; i++) {
    int64_t eindex, mmp_index;
    eindex = id_to_index(evolved, now.halos[i].descid);
    if (eindex < 0) continue;
    if (evolved.halos[eindex].mmp_id == now.halos[i].id) continue;
    mmp_index = id_to_index(now, evolved.halos[eindex].mmp_id);
    if (mmp_index < 0) continue;
    if (now.halos[mmp_index].mvir * MAJOR_MERGER < now.halos[i].mvir)
      evolved.halos[eindex].flags |= MAJOR_MERGER_FLAG;
  }
}

void calculate_tracking_stats(float a1, float a2) {
  int64_t i, eindex;
  //#pragma omp parallel for private(i,eindex) schedule(static,1000)
  for (i=0; i<now.num_halos; i++) {
    eindex = id_to_index(evolved, now.halos[i].descid);
    if (eindex < 0) continue;
    if (really_beyond_mmp_ratio(now.halos[i], evolved.halos[eindex])) continue;
    now.halos[i].tracked = evolved.halos[eindex].tracked + 1;

    if (!(evolved.halos[eindex].flags & MAJOR_MERGER_FLAG) && 
	now.halos[i].id == evolved.halos[eindex].mmp_id) {
      now.halos[i].tracked_single_mmp = evolved.halos[eindex].tracked_single_mmp;
      if (!now.halos[i].phantom) now.halos[i].tracked_single_mmp++;
    }
    now.halos[i].num_mmp_phantoms = evolved.halos[eindex].num_mmp_phantoms;
    if (now.halos[i].phantom > 0) now.halos[i].num_mmp_phantoms++;
    if (((1+now.halos[i].tracked)*MAX_PHANTOM_FRACTION < 
         (now.halos[i].num_mmp_phantoms-MAX_PHANTOM)) && 
        (now.halos[i].tracked >= MIN_TIMESTEPS_TRACKED)) {
      now.halos[i].flags |= DEAD_HALO_FLAG;
      log_too_many_phantoms(a1, a2, &(now.halos[i]));
    }
  }
}

void create_phantom_progenitors(float a1, float a2)
{
  int64_t i, j, last_id = 0, index, eindex;
  struct tree_halo tmp_halo;
  float *incoming_mass = NULL;
  incoming_mass = check_realloc(incoming_mass, sizeof(float)*evolved.num_halos,
				"Allocating incoming mass checks.");
  memset(incoming_mass, 0, sizeof(float)*evolved.num_halos);
  //#pragma omp parallel for private(i,eindex)
  for (i=0; i<now.num_halos; i++) {
    eindex = id_to_index(evolved, now.halos[i].descid);
    if (eindex < 0) continue;
    //#pragma omp atomic
    incoming_mass[eindex]+=now.halos[i].mvir;
  }

  for (j=0; j<evolved.num_halos; j++) {
    index = id_to_index(now, evolved.halos[j].mmp_id);
    if (index >= 0) { //Create phantoms for halos with apparently incorrectly-identified progenitors.
      if (evolved.halos[j].mvir*MIN_MMP_MASS_RATIO < incoming_mass[j])
	continue;
      if (!beyond_mmp_ratio(now.halos[index], evolved.halos[j]))
	continue;
    }
    tmp_halo = evolved.halos[j];
    tmp_halo.phantom = evolved.halos[j].phantom + 1;
    if (tmp_halo.phantom > MAX_PHANTOM ||
	(tmp_halo.phantom > MAX_PHANTOM_SMALL && 
	 tmp_halo.np < SMALL_PARTICLE_LIMIT)) {
      //Flag halo for removal
      tmp_halo.flags |= DEAD_HALO_FLAG;
      log_dead_halo(a1, &tmp_halo);
    }
    last_id = gen_new_phantom_id(&now, last_id);
    tmp_halo.id = last_id;
    tmp_halo.descid = evolved.halos[j].id;
    tmp_halo.pid = tmp_halo.upid = tmp_halo.tidal_id = tmp_halo.mmp_id = -1;
    evolved.halos[j].mmp_id = tmp_halo.id;
    add_halo(&now, tmp_halo);
    log_phantom_halo(a1, a2, &tmp_halo, &(evolved.halos[j]));
  }
  build_id_conv_list(&now);

  free(incoming_mass);
}

void mark_dead_halos(float a) {
  int64_t j;
  //Mark all remaining halos without descendants as dead.
#pragma omp parallel for private(j) schedule(static,1000)
  for (j=0; j<now.num_halos; j++) {
    int64_t index = id_to_index(evolved, now.halos[j].descid);
    if (index < 0) {
      now.halos[j].flags |= DEAD_HALO_FLAG;
      log_dead_halo(a, &(now.halos[j]));
    }
  }
}

void fix_halos_without_descendants(float a1, float a2) {
  int64_t j,k;
  
  //Assign halos w/o descendants based on tidal forces
  //Have to do this recursively in case nested mergers occur 
  for (k=0; k<RECURSION_LIMIT; k++) {
    //#pragma omp parallel for private(j)
    for (j=0; j<now.num_halos; j++) {
      int64_t index, pindex, eindex;
      if (now.halos[j].descid >= 0) continue;
      if (now.halos[j].tidal_id < 0) continue;
      index = id_to_index(now, now.halos[j].tidal_id);
      if (index < 0 || now.halos[j].tidal_force < TIDAL_FORCE_LIMIT) {
	now.halos[j].tidal_id = -1;
	continue;
      }
      now.halos[j].descid = now.halos[index].descid;
      if (now.halos[j].descid<0) continue;
      eindex = id_to_index(evolved, now.halos[index].descid);
      if (eindex < 0) continue;
      if (evolved.halos[eindex].mvir < (now.halos[j].mvir + now.halos[index].mvir)*MIN_MMP_MASS_RATIO) {
	fprintf(logfile, "Tidal mass conservation impossible (%"PRId64" %e + %"PRId64" %e != %"PRId64" %e)\n", now.halos[j].id, now.halos[j].mvir,
		now.halos[index].id, now.halos[index].mvir, evolved.halos[eindex].id, evolved.halos[eindex].mvir);
	pindex = id_to_index(evolved, evolved.halos[eindex].upid);
	if (pindex >= 0 && !(evolved.halos[pindex].flags & DEAD_HALO_FLAG)) {
	  now.halos[j].descid = evolved.halos[pindex].id;
	  fprintf(logfile, "Tidal mass stripping exception: %"PRId64" -> %"PRId64"\n", now.halos[j].id, evolved.halos[pindex].id);
	} else {
	  pindex = id_to_index(now, now.halos[index].upid);
	  if (pindex < 0) pindex = id_to_index(now, now.halos[j].upid);
	  if (pindex >= 0 && now.halos[pindex].descid>=0 && !(now.halos[pindex].flags & DEAD_HALO_FLAG)) {
	    now.halos[j].descid = now.halos[pindex].descid;
	    fprintf(logfile, "Tidal mass stripping exception (2): %"PRId64" -> %"PRId64"\n", now.halos[j].id, now.halos[pindex].descid);
	  }
	}
      }
      log_tidal_repair(a1, a2, &(now.halos[j]), &(now.halos[index]));
    }
  }
}


void fix_halos_without_progenitors(float a1, float a2, int64_t stage) {
  int64_t i, j, k, index, num_matches;
  float d, min_metric;
  struct fast3tree_results *nearest;
  struct tree_halo *h2;
  struct tree_halo tmp_halo;
  float sigma_x, sigma_v;
  float vmax_tol, vmax_avg;
  float metric_limit = (stage > 1) ? 2*METRIC_LIMIT : METRIC_LIMIT;
  float dist_limit = 0;
  float max_dist = (float) box_size / 2.01;

  //Build tree of halos without descendants
  build_nodesc_halo_tree(now.halos, now.num_halos);
  nearest = fast3tree_results_init();

  //Find halos without progenitors:
  for (i=0; i<evolved.num_halos; i++) {
    //j = (stage == 3) ? process_order[i] : i;
    j = process_order[i];
    if (evolved.halos[j].mmp_id>=0) continue;
    if ((stage == 2) && (evolved.halos[j].mvir < MASS_RES_OK)) continue;

    sigma_x_v_vmax(&(evolved.halos[j]), &sigma_x, &sigma_v, &vmax_avg, &vmax_tol, a1);
    dist_limit = (stage >= 3) ? 
      ((LAST_DITCH_SEARCH_LIMIT*1e-3)*evolved.halos[j].rvir) 
      : (sigma_x*metric_limit);
    if (dist_limit > max_dist) dist_limit = max_dist;

    IF_PERIODIC
      fast3tree_find_sphere_periodic(halo_tree, nearest, evolved.halos[j].pos,
				     dist_limit);
    else
      fast3tree_find_sphere(halo_tree, nearest, evolved.halos[j].pos,
			    dist_limit);
    
    min_metric = metric_limit;
    num_matches = 0;
    for (k=0; k<nearest->num_points; k++) {
      h2 = nearest->points[k]->h;
      if (h2->descid >= 0) continue; //Ignore halos already with descendants
      if ((stage <= 2) && beyond_mmp_ratio(*h2, evolved.halos[j])) continue;
      if ((stage == 2) && (h2->mvir < MASS_RES_OK)) continue;
      d = (stage >= 3) ? last_ditch_metric(&(evolved.halos[j]), h2)
	: metric(h2, &(evolved.halos[j]), sigma_x, sigma_v, vmax_avg, vmax_tol);

      if ((stage == 3) && (d > LAST_DITCH_VMAX_RATIO_1)) {
	if (evolved.halos[j].phantom < MAX_PHANTOM_SMALL)
	  d = UNPHYSICAL;
      }

      if (d < min_metric) {
	min_metric = d;
	evolved.halos[j].mmp_id = h2->id;
      }
      if (d < metric_limit) {
	num_matches++;
	if (num_matches == 1)
	  fprintf(logfile, "Progenitor candidates for %f %"PRId64" %e %f: ", a2,
		  evolved.halos[j].id, evolved.halos[j].mvir, evolved.halos[j].vmax);
	if (stage < 3)
	  fprintf(logfile, "%"PRId64" %e %f (%f) (stage %"PRId64") ", h2->id, h2->mvir, h2->vmax, d, stage);
      }
    }

    if (num_matches) fprintf(logfile, "(%"PRId64" match(es))\n", num_matches);

    if (evolved.halos[j].mmp_id<0) continue;
    index = id_to_index(now, evolved.halos[j].mmp_id);
    assert(index >= 0);
    now.halos[index].descid = evolved.halos[j].id;
    if (stage > 1) now.halos[index].flags |= SUSPICIOUS_LINK_FLAG;
    if (stage > 2) now.halos[index].flags |= UNPHYSICAL_LINK_FLAG;
    log_grav_repair(a1, a2, min_metric, &(now.halos[index]), &(evolved.halos[j]));
  }
  fast3tree_results_free(nearest);

  //Copy actual halo properties for evolved halos with identified MMPs:
#pragma omp parallel for private(j,index,tmp_halo) schedule(static,1000)
  for (j=0; j<evolved.num_halos; j++) {
    if (evolved.halos[j].mmp_id < 0) continue;
    index = id_to_index(now, evolved.halos[j].mmp_id);
    assert(index >= 0);
    if (beyond_mmp_ratio(now.halos[index], evolved.halos[j])) continue;
    tmp_halo = evolved.halos[j];
    evolved.halos[j] = now.halos[index];
    evolved.halos[j].id = tmp_halo.id;
    evolved.halos[j].mvir = tmp_halo.mvir;
    evolved.halos[j].vmax = tmp_halo.vmax;
    evolved.halos[j].pid = tmp_halo.pid;
    evolved.halos[j].upid = tmp_halo.upid;
    evolved.halos[j].descid = tmp_halo.descid;
    evolved.halos[j].mmp_id = tmp_halo.mmp_id;
    evolved.halos[j].flags = tmp_halo.flags;
    evolved.halos[j].phantom = tmp_halo.phantom;
    evolved.halos[j].tracked = tmp_halo.tracked;
    evolved.halos[j].tracked_single_mmp = tmp_halo.tracked_single_mmp;
    evolved.halos[j].num_mmp_phantoms = tmp_halo.num_mmp_phantoms;
    evolved.halos[j].spin = tmp_halo.spin;
    memcpy(evolved.halos[j].J, tmp_halo.J, sizeof(float)*3);
  }
}

void mark_unphysical_halos(float a) {
  int64_t j;
#pragma omp parallel for private(j) schedule(static,1000)
  for (j=0; j<now.num_halos; j++) {
    if ((now.halos[j].mvir <=0) || (now.halos[j].vmax <= 0)) {
      fprintf(logfile, "Unphysical halo: marking for deletion. (%f %"PRId64" %e %f)\n",
	      a, now.halos[j].id, now.halos[j].mvir, now.halos[j].vmax);
      now.halos[j].descid = -1;
      now.halos[j].tidal_force = 0;
      now.halos[j].tidal_id = -1;
      now.halos[j].flags |= DEAD_HALO_FLAG;
      continue;
    }
  } 
}

void break_spurious_links(float a1, float a2) {
  int64_t j;
  //Break spurious links in the merger tree
  //#pragma omp parallel for private(j)
  for (j=0; j<now.num_halos; j++) {
    float sigma_x, sigma_v, sigma_vmax, d, vmax_avg;
    int64_t eindex = id_to_index(evolved, now.halos[j].descid);
    if (eindex < 0) {
      if (now.halos[j].descid >= 0) log_desc_not_found(a1, a2, &(now.halos[j]));
      now.halos[j].descid = -1;
      continue;
    }

    //Check to see whether we're the MMP or not
    if (evolved.halos[eindex].mmp_id != now.halos[j].id) {
      //If we're not, check the tidal force estimate
      if (now.halos[j].tidal_force < TIDAL_FORCE_LIMIT) {
	now.halos[j].descid = -1;
	log_not_enough_tidal(a1, a2, &(now.halos[j]), &(evolved.halos[eindex]));
      }
      continue;
    }

    //Otherwise, reject link if the evolved halo does not match the MMP in the tree
    sigma_x_v_vmax(&(evolved.halos[eindex]), &sigma_x, &sigma_v, &vmax_avg, &sigma_vmax, a1);
    d = metric(&(now.halos[j]), &(evolved.halos[eindex]), sigma_x, sigma_v, vmax_avg, sigma_vmax);
    if (d > METRIC_BREAK_LIMIT) {
      now.halos[j].descid = -1;
      evolved.halos[eindex].mmp_id = -1;
      log_grav_inconsistency(a1, a2, d, &(now.halos[j]), &(evolved.halos[eindex]));
    }
  }
}


void break_non_mmp_links(float a1, float a2) {
  int64_t j, eindex;

  for (j=0; j<now.num_halos; j++) {
    eindex = id_to_index(evolved, now.halos[j].descid);
    if (eindex < 0) continue;

    if (now.halos[j].id != evolved.halos[eindex].mmp_id) {
      now.halos[j].descid = -1; //Break non-mmp links.
      continue;
    }
    
    if (beyond_mmp_ratio(now.halos[j], evolved.halos[eindex])) {
      log_spurious_mmp_link_broken(a1, a2, &(now.halos[j]), &(evolved.halos[eindex]));
      now.halos[j].descid = -1;
      evolved.halos[eindex].mmp_id = -1;
    }
  }
}

void build_mmp_list(int64_t stage, float a1, float a2)
{
  int64_t j, eindex, index;
  //Build up the most-massive progenitor list
  //This simultaneously tags halos w/o progenitors
  for (j=0; j<evolved.num_halos; j++) evolved.halos[j].mmp_id = -1;

  for (j=0; j<now.num_halos; j++) {
    if (now.halos[j].flags & DEAD_HALO_FLAG) {
      now.halos[j].flags |= MMP_FLAG;
      continue;
    }
 
    eindex = id_to_index(evolved, now.halos[j].descid);
    if (stage == 1 && eindex < 0) { //Descendant halo no longer exists
      if (now.halos[j].descid >= 0) log_desc_not_found(a1, a2, &(now.halos[j]));
      else log_no_desc(a1, &(now.halos[j]));
      now.halos[j].descid = -1;
      continue;
    } else {
      assert(eindex>=0); //All links should now be fixed.
      //Not an actual link if real descendant was 
      //deleted and replaced by a phantom
      if ((stage == 1) && (evolved.halos[eindex].phantom>0)) { 
	log_desc_not_found(a1, a2, &(now.halos[j]));
	now.halos[j].descid = -1;
	continue;
      }
    }

    index = id_to_index(now, evolved.halos[eindex].mmp_id);
    if (index < 0) { //If no MMP already identified.
      evolved.halos[eindex].mmp_id = now.halos[j].id;
      if (stage != 1) now.halos[j].flags |= MMP_FLAG;
      continue;
    }

    //Otherwise, check to see if we're more massive than the previous MMP.
    if (now.halos[j].mvir >= now.halos[index].mvir) {
      evolved.halos[eindex].mmp_id = now.halos[j].id;
      if (stage != 1) {
	now.halos[j].flags |= MMP_FLAG;
	now.halos[index].flags -= (now.halos[index].flags & MMP_FLAG); //Safer than &= ~MMP_FLAG
      }
    }
  }
}


void initialize_last_timestep(int64_t last_output_num, float scale)
{
  char buffer[1024];
  int64_t i;
  timed_output(_GC, "Loading snapshot %"PRId64"...", last_output_num);
  snprintf(buffer, 1024, "%s/out_%"PRId64".list%s", INBASE, last_output_num,
	   ((!strcasecmp(INPUT_FORMAT, "BINARY")) ? ".bin" : ""));
  gc_load_halos(buffer, &now, scale, 1);
  output_order = check_realloc(output_order, now.num_halos*sizeof(int64_t), "halo output order");

  for (i=0; i<now.num_halos; i++) {
    now.halos[i].flags |= MMP_FLAG;
    output_order[i] = now.halos[i].id;
  }

  copy_halos_to_evolved(last_output_num, 0);
}

void regen_output_order(int64_t output_num, float scale) {
  char buffer[1024];
  int64_t i;
  timed_output(_GC, "Reloading snapshot %"PRId64"...", output_num);
  snprintf(buffer, 1024, "%s/consistent_%"PRId64".list", OUTBASE, output_num);
  gc_load_halos(buffer, &evolved, scale, 0);
  output_order = check_realloc(output_order, evolved.num_halos*sizeof(int64_t), "halo output order");
  for (i=0; i<evolved.num_halos; i++)
    output_order[i] = evolved.halos[i].id;
}

int desc_order_compare(const void *a, const void *b)
{
  int64_t c = id_to_index(now, *((int64_t *)a));
  int64_t d = id_to_index(now, *((int64_t *)b));
  assert(c>=0 && d>=0);
  int64_t dc = now.halos[c].descid;
  int64_t dd = now.halos[d].descid;
  float mass_compare = now.halos[c].mvir - now.halos[d].mvir;
  if (dc<evolved.min_id) {
    if (dd<evolved.min_id) {
      if (mass_compare > 0) return 1; //Sort in increasing order of mass
      if (mass_compare < 0) return -1;
      return 0;
    }
    return 1;
  }
  if (dd<evolved.min_id) return -1;
  assert((dd<=evolved.max_id) && (dc<=evolved.max_id));
  int64_t dco = output_order_inv[id_to_index(evolved, dc)];
  int64_t ddo = output_order_inv[id_to_index(evolved, dd)];
  int64_t id_compare = dco - ddo;
  if (id_compare > 0) return 1;
  if (id_compare < 0) return -1;
  if (mass_compare > 0) return -1; //Sort in decreasing order of mass
  if (mass_compare < 0) return 1;
  return 0;
}

void calculate_new_output_order(void)
{
  int64_t i, eindex;
  output_order_inv=check_realloc(output_order_inv, sizeof(int64_t)*evolved.num_halos, "output order inversion");

  for (i=0; i<evolved.num_halos; i++) {
    eindex = id_to_index(evolved, output_order[i]);
    assert(eindex>=0);
    output_order_inv[eindex] = i;
  }

  output_order = check_realloc(output_order, sizeof(int64_t)*now.num_halos, "new output order");
  for (i=0; i<now.num_halos; i++) output_order[i] = now.halos[i].id;

  qsort(output_order, now.num_halos, sizeof(int64_t), desc_order_compare);
}

void fork_and_print_halos(char *buffer, struct tree_halo *h, int64_t num_halos, int64_t gzip_halos) {
  struct tree_halo *my_halos = NULL;
  pid_t pid;
  
  my_halos = check_realloc(NULL, sizeof(struct tree_halo)*num_halos,
			   "Allocating output halo memory.");
  memcpy(my_halos, h, sizeof(struct tree_halo)*num_halos);

  pid = fork();
  if (pid) free(my_halos);
  if (pid < 0) common_print_halos(buffer, h, num_halos, gzip_halos);
  if (pid) return;

  clear_halo_stash(&now);
  clear_halo_stash(&evolved);
  if (halo_tree) fast3tree_rebuild(halo_tree, 0, NULL);
  if (halo_positions) free(halo_positions);
  free_tidal_trees();
  common_print_halos(buffer, my_halos, num_halos, gzip_halos);
  exit(0);
}


void copy_halos_to_evolved(int64_t output_num, int64_t prune_dead_halos)
{
  int64_t i, index, num_good_halos=0, num_dead_halos=0;
  struct tree_halo *th;
  char buffer[1024];
  pid_t pid;

  evolved.num_halos = now.num_halos;
  evolved.halos = check_realloc(evolved.halos,
	sizeof(struct tree_halo)*now.num_halos, "Allocating halos.");

  for (i=0; i<now.num_halos; i++) {
    index = id_to_index(now, output_order[i]);
    th = &(now.halos[index]);

    if (prune_dead_halos && (th->flags & DEAD_HALO_FLAG)) {
      num_dead_halos++;
      memcpy(&(evolved.halos[evolved.num_halos - num_dead_halos]), th,
	     sizeof(struct tree_halo));
    } else {
      th->flags -= (th->flags & DEAD_HALO_FLAG);
      memcpy(&(evolved.halos[num_good_halos]), th, sizeof(struct tree_halo));
      if (num_good_halos < i)
	output_order[num_good_halos] = output_order[i];
      num_good_halos++;
    }
  }
  if (logfile) {
    fprintf(logfile, "#Good halos: %"PRId64"\n#Dead halos: %"PRId64"\n", num_good_halos, num_dead_halos);
    fclose(logfile);

    pid = fork();
    if (pid <= 0) {
      snprintf(buffer, 1024, "%s/logfile_%"PRId64".list", OUTBASE, output_num);
      gzip_file(buffer);
      if (!pid) exit(0);
    }
  }

  output_order = check_realloc(output_order, sizeof(int64_t)*num_good_halos, "output order");
  clear_halo_stash(&now);

  snprintf(buffer, 1024, "%s/consistent_%"PRId64".list", OUTBASE, output_num);
  fork_and_print_halos(buffer, evolved.halos, num_good_halos, 0);
  snprintf(buffer, 1024, "%s/dead_%"PRId64".list", OUTBASE, output_num);
  fork_and_print_halos(buffer, evolved.halos + num_good_halos, num_dead_halos, 0);

  evolved.num_halos = num_good_halos;
  evolved.halos = check_realloc(evolved.halos,
	sizeof(struct tree_halo)*evolved.num_halos, "Allocating halos.");
  build_id_conv_list(&evolved);
}

void clear_halo_stash(struct halo_stash *h) {
  h->halos = (struct tree_halo *)realloc(h->halos, 0);
  h->id_conv = (int64_t *)realloc(h->id_conv, 0);
  h->max_id = h->min_id = h->num_halos = 0;
}

int sort_evolved_halo_list_by_decreasing_vmax(const void *a, const void *b) {
  struct halo_stash *h = &evolved;
  struct tree_halo *c = h->halos + *((int64_t *)a);
  struct tree_halo *d = h->halos + *((int64_t *)b);
  if (c->vmax > d->vmax) return -1;
  if (c->vmax < d->vmax) return 1;
  return 0;
}

void sort_process_list_by_vmax_evolved() {
  int64_t i;
  struct halo_stash *h = &evolved;
  process_order = check_realloc(process_order, sizeof(int64_t)*h->num_halos, "process order");
  for (i=0; i<h->num_halos; i++) process_order[i] = i;
  qsort(process_order, h->num_halos, sizeof(int64_t), sort_evolved_halo_list_by_decreasing_vmax);
}

void gc_load_halos(char *filename, struct halo_stash *h, float scale, int64_t orig_data)
{ 
  int64_t i;
  box_size = 0;
  clear_halo_stash(h);
  load_halos(filename, h, scale, 0);
  build_id_conv_list(h);
  if (orig_data && (FIX_ROCKSTAR_SPINS > -1)) {
    assert(FIX_ROCKSTAR_SPINS < EXTRA_PARAMS);
    for (i=0; i<h->num_halos; i++) {
      h->halos[i].extra_params[FIX_ROCKSTAR_SPINS] *= 2.0;
      float t_u = h->halos[i].extra_params[FIX_ROCKSTAR_SPINS];
      if (fabs(2.0-t_u) != 0)
	h->halos[i].spin *= sqrt(fabs(1.0-t_u)/fabs(2.0-t_u));
    }
  }
}

/* Max range is the maximum distance over which the gravitational effects of
   the largest halo will matter. */
/* Exclude halos already with descendants. */
void build_nodesc_halo_tree(struct tree_halo *h, int64_t num_h) {
  int64_t i,j,num_pos = 0;

  halo_positions = (struct halo_pos *)check_realloc(halo_positions,
      sizeof(struct halo_pos)*num_h, "Allocating halo positions.");

  for (i=0; i<num_h; i++) {
    if (h[i].descid >= 0) continue;
    halo_positions[num_pos].h = &(h[i]);
    for (j=0; j<3; j++)
      halo_positions[num_pos].pos[j] = h[i].pos[j];
    num_pos++;
  }

  halo_positions = (struct halo_pos *)check_realloc(halo_positions,
      sizeof(struct halo_pos)*num_pos, "Allocating halo positions.");

  fast3tree_rebuild(halo_tree, num_pos, halo_positions);
  IF_PERIODIC _fast3tree_set_minmax(halo_tree, 0, box_size);
}


void calc_metric_stats(int64_t output_num, float a1, float a2) {
  int64_t j, eindex, index;
  for (j=0; j<now.num_halos; j++) {
    if (now.halos[j].flags & DEAD_HALO_FLAG)
      continue;
 
    eindex = id_to_index(evolved, now.halos[j].descid);
    if (eindex < 0) //Descendant halo no longer exists
      continue;
    if (evolved.halos[eindex].mmp_id == now.halos[j].id) {
      index = id_to_index(now, now.halos[j].pid);
      if (index < 0)
	build_metric_stats(&(now.halos[j]), NULL, &(evolved.halos[eindex]));
      else
	build_metric_stats(&(now.halos[j]), &(now.halos[index]), 
			   &(evolved.halos[eindex]));
    }
  }
  finish_metric_stats(output_num, a1, a2);
}


float last_ditch_metric(struct tree_halo *h1, struct tree_halo *h2)
{
    if (!(h1->mvir > 0) || !(h2->mvir > 0)
      || !(h1->vmax > 0) || !(h2->vmax > 0)) return UNPHYSICAL;
    float vmax_ratio = (h1->vmax > h2->vmax) ? (h1->vmax / h2->vmax) :
      (h2->vmax / h1->vmax);
    if (vmax_ratio > LAST_DITCH_VMAX_RATIO_2) return UNPHYSICAL;
    return vmax_ratio;
}

float metric(struct tree_halo *h1, struct tree_halo *h2, float sigma_x, float sigma_v, float vmax_avg, float vmax_tol)
{
  //dt is in Myr
  int64_t i;
  float dx, dr=0;
  float dvx, dv=0;
  float dvmax = 0;
  
  if (!(h1->mvir > 0) || !(h2->mvir > 0)
      || !(h1->vmax > 0) || !(h2->vmax > 0)) return UNPHYSICAL;
  dvmax = fabsf((log10f(h1->vmax/h2->vmax)-vmax_avg) / vmax_tol);
  
  for (i=0; i<3; i++) {
    dx = h1->pos[i] - h2->pos[i];
    IF_PERIODIC {
      if (dx>box_size/2.0) dx -= box_size;
      if (dx<-box_size/2.0) dx += box_size;
    }
    dr+=dx*dx;

    dvx = h1->vel[i] - h2->vel[i];
    dv+=dvx*dvx;
  }

  assert((sigma_x > 0) && (sigma_v > 0));
  return (sqrtf(dr/(2.0*sigma_x*sigma_x) + dv/(2.0*sigma_v*sigma_v)
	       + dvmax*dvmax/2.0));
}

void print_tidal_forces(int64_t output_num, float a, struct halo_stash *h, struct halo_stash *evolved) {
  char buffer[1024];
  FILE *mmp_out, *nmmp_out, *output;
  int64_t mmp, j, eindex, tidal_desc_id, nindex;
  float tidal_mvir;
  snprintf(buffer, 1024, "%s/tidal_forces_mmp_%"PRId64".list", OUTBASE, output_num);
  mmp_out = check_fopen(buffer, "w");
  fprintf(mmp_out, "#ID DescID TidalId TidalDescId TidalForce Rvir Mvir TidalMvir Vmax\n");
  fprintf(mmp_out, "#Scale factor: %f\n", a);
  snprintf(buffer, 1024, "%s/tidal_forces_nmmp_%"PRId64".list", OUTBASE, output_num);
  nmmp_out = check_fopen(buffer, "w");
  fprintf(nmmp_out, "#ID DescID TidalId TidalDescId TidalForce Rvir Mvir TidalMvir Vmax\n");
  fprintf(nmmp_out, "#Scale factor: %f\n", a);
  for (j=0; j<h->num_halos; j++) {
    if (h->halos[j].flags & DEAD_HALO_FLAG)
      continue;
 
    eindex = id_to_index(*evolved, h->halos[j].descid);
    if (eindex < 0) //Descendant halo no longer exists
      continue;
    mmp = (evolved->halos[eindex].mmp_id == h->halos[j].id) ? 1 : 0;
    //if (mmp && really_beyond_mmp_ratio(h->halos[j], evolved->halos[eindex])) mmp = 0;
    output = (mmp) ? mmp_out : nmmp_out;

    nindex = id_to_index(*h, h->halos[j].tidal_id);
    if (nindex>-1) {
      tidal_desc_id = h->halos[nindex].descid;
      tidal_mvir = h->halos[nindex].mvir;
    }
    else { tidal_desc_id = -1; tidal_mvir = 0; }

    fprintf(output, "%"PRId64" %"PRId64" %"PRId64" %"PRId64" %.3e %f %.3e %.3e %f\n",
	    h->halos[j].id, h->halos[j].descid, h->halos[j].tidal_id,
	    tidal_desc_id, h->halos[j].tidal_force, h->halos[j].rvir, 
	    h->halos[j].mvir, tidal_mvir, h->halos[j].vmax);
  }
  fclose(mmp_out);
  fclose(nmmp_out);
}
