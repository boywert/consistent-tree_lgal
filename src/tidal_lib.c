#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "distance.h"
#include "universe_time.h"
#include "masses.h"
#include "halo_evolve_lib.h"
#include "tree_halo.h"
#include "gravitational_consistency.h"
#include "check_syscalls.h"
#include "spinlock.h"

#define FAST3TREE_PREFIX TIDAL_LIB
#define FAST3TREE_TYPE struct halo_pos 
#include "fast3tree.c"

struct halo_pos *hpos = NULL;
struct fast3tree *tidal_tree = NULL;
extern float box_size; /* Automatically set; in Mpc/h */
extern double h0; /* h0 */
int64_t tidal_extra_range = 0;


void free_tidal_trees(void) {
  if (tidal_tree) fast3tree_rebuild(tidal_tree, 0, NULL);
  if (hpos) free(hpos);
}

/* Assumes mvir in Msun and dt in Myr. */
inline float halo_tidal_range(float mvir, float av_a) {
  /* TIDAL_FORCE_LIMIT = Gc*m/r^3 (km/s/Myr per comoving Mpc) */
  /* r^3 = Gc*m/TIDAL_FORCE_LIMIT / a^2 */
  /* Mutiply by 2 to be safe.*/
  float tforce_limit = TIDAL_FORCE_LIMIT;
  if (tidal_extra_range && TIDAL_FORCE_LIMIT > 0.001) {
    tforce_limit = 0.01 * TIDAL_FORCE_LIMIT;
    if (tforce_limit < 0.001) tforce_limit = 0.001;
  }
  return (2*cbrtf(fabs(mvir*Gc/(tforce_limit*av_a*av_a))));
}

void calc_tidal_forces(struct halo_stash *h, double a1, double a2)
{
  int64_t i;
  float dt = (scale_to_years(a2) - scale_to_years(a1))/1.0e6; //In Myr
  float av_a = (a1+a2)/2.0;
  float acorr = 1; //pow(av_a, 0.7);
  float conv_const = av_a*av_a / h0 / h0 / h0;
  float vel_dt = dt * 1.02268944e-6 * h0 / av_a; //1 km/s to comoving Mpc/Myr/h
  struct tree_halo *halos = h->halos;
  float max_dist = box_size / 2.0;
  spinlock_t * volatile locks;
  ALLOC_SPINLOCKS(locks, h->num_halos);

  gen_ff_cache();

  if (!tidal_tree) tidal_tree = fast3tree_init(0, NULL);

  hpos = check_realloc(hpos, sizeof(struct halo_pos)*h->num_halos, "Allocating tidal halo positions");
  vel_dt /= 2.0;
#pragma omp parallel for private(i)
  for (i=0; i<h->num_halos; i++) {
    halos[i].tidal_force = 0;
    halos[i].tidal_id = -1;
    halos[i].pid = halos[i].upid = -1;
    halos[i].pid_vmax = halos[i].upid_vmax = 0;
    hpos[i].h = &(halos[i]);
    for (int64_t j=0; j<3; j++) {
      hpos[i].pos[j] = halos[i].pos[j] + vel_dt*halos[i].vel[j];
      IF_PERIODIC {
	if (hpos[i].pos[j] > box_size)
	  hpos[i].pos[j] -= box_size;
	else if (hpos[i].pos[j] < 0)
	  hpos[i].pos[j] += box_size;
      }
    }
  }
  fast3tree_rebuild(tidal_tree, h->num_halos, hpos);
  IF_PERIODIC _fast3tree_set_minmax(tidal_tree, 0, box_size);
  
#pragma omp parallel default(none) private(i) shared(hpos,h,h0,av_a,max_dist,tidal_tree,locks,MAJOR_MERGER,acorr,box_size,conv_const,halos)
  {
    //Calculate accelerations
    struct fast3tree_results *nearest = fast3tree_results_init();
#pragma omp for
    for (i=0; i<h->num_halos; i++) {
      struct tree_halo *h1 = hpos[i].h;
      if (h1->flags & MERGER_FLUCTUATION_FLAG) continue;
      float range = halo_tidal_range(h1->mvir/h0, av_a)*h0; //In comoving Mpc/h
      float rvir = h1->rvir / 1.0e3; //In comoving Mpc/h
      float nearest_range = rvir;
      if (range > rvir) nearest_range = range;
      IF_PERIODIC {
	if (nearest_range > max_dist) nearest_range = max_dist;
	fast3tree_find_sphere_periodic(tidal_tree, nearest, hpos[i].pos, nearest_range);
      } else {
	fast3tree_find_sphere(tidal_tree, nearest, hpos[i].pos, nearest_range);
      }
      
      float inv_rs = 1000.0/h1->rs; //In comoving h/Mpc
      float mass_factor = calculate_mass_factor(h1->mvir/h0, h1->rvir, h1->rs);
      for (int64_t j=0; j<nearest->num_points; j++) {
	//Calculate tidal forces
	struct tree_halo *h2 = nearest->points[j]->h;
	
#ifndef NO_PERIODIC
#define DIST(a,b) float a = hpos[i].pos[b] - nearest->points[j]->pos[b]; \
	if (a > max_dist) a-=box_size;					\
	else if (a < -max_dist) a+=box_size;
#else
#define DIST(a,b) float a = hpos[i].pos[b] - nearest->points[j]->pos[b];
#endif
	DIST(dx,0);
	DIST(dy,1);
	DIST(dz,2);
#undef DIST
	float r = sqrtf(dx*dx + dy*dy + dz*dz); // in comoving Mpc/h
      
	//Including acorr slightly improves velocity results
	// as a function of redshift.
	float m = mass_factor*ff_cached(r*inv_rs*acorr); //In Msun
	
	float r3 = r * r * r * conv_const; //in (real Mpc)^3
	float acc = r3 ? Gc*m/r3 : 0; //In km/s / Myr / (real Mpc)
	
	int64_t index = h2-halos;
	assert(index > -1 && index < h->num_halos);
	acquire_spinlock(locks, index);
	if (r < rvir && h1->vmax > h2->vmax) {
	  if (h2->pid < 0) {
	    h2->pid = h2->upid = h1->id;
	    h2->pid_vmax = h2->upid_vmax = h1->vmax;
	  }
	  else {
	    if (h2->pid_vmax > h1->vmax) {
	      h2->pid = h1->id;
	      h2->pid_vmax = h1->vmax;
	    }
	    if (h2->upid_vmax < h1->vmax) {
	      h2->upid = h1->id;
	      h2->upid_vmax = h1->vmax;
	    }
	  }
	}
	
	if ((acc > h2->tidal_force) && ((h2->mvir*MAJOR_MERGER) < h1->mvir)
	    && (!(h2->phantom))) {
	  h2->tidal_force = acc;
	  h2->tidal_id = h1->id;
	}
	release_spinlock(locks, index);
      }
    }

    fast3tree_results_free(nearest);
  }
  free(locks);
}
