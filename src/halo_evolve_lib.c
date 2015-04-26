#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include "distance.h"
#include "tree_halo.h"
#include "masses.h"
#include "universe_time.h"
#include "halo_evolve_lib.h"
#include "stringparse.h"
#include "gravitational_consistency.h"
#include "halo_io.h"

#define FAST3TREE_TYPE struct tree_halo 
#define FAST3TREE_PREFIX HEVOLVE
#include "fast3tree.c"

struct fast3tree *ev_halo_tree = NULL;
extern float box_size; /* Automatically set; in Mpc/h */
extern float max_mvir; /* Automatically set; in Msun */
extern double h0; /* Hubble constant */


void evolve_halos(float a1, float a2, struct halo_stash *h) {
  float da, a;
  int64_t i;
  struct tree_halo *h1;

  //#pragma omp parallel for private(i,h1)
  for (i=0; i<h->num_halos; i++) {
    h1 = &(h->halos[i]);
    h1->mass_factor = calculate_mass_factor(h1->mvir /h0, h1->rvir, h1->rs);
  }

  if (!ev_halo_tree) ev_halo_tree = fast3tree_init(0, h->halos);
  da = (a2-a1)/((double)NUM_TIMESTEPS);
  do_timestep(h, 0, a1, a1+da, 1); //First step
  for (i=1; i<NUM_TIMESTEPS; i++) {
    a = a1+i*da;
    do_timestep(h, a-da, a, a+da, 0);
  }
  do_timestep(h, a, a+da, 0, -1); //Last step
}

void correct_mass_factors(struct halo_stash *h) {
#pragma omp parallel default(none) shared(h,box_size,ev_halo_tree)
  {
    int64_t i,j;
    float max_dist = box_size / 2.0;
    float sat_mass, range;
    struct fast3tree_results *nearest = fast3tree_results_init();
#pragma omp for
    for (i=0; i<h->num_halos; i++) {
      range = h->halos[i].rvir/1e3;
      IF_PERIODIC {
	if (max_dist < range) range = max_dist*0.99;
	assert(fast3tree_find_sphere_periodic(ev_halo_tree, nearest, 
					    h->halos[i].pos, range));
      }
      else
	fast3tree_find_sphere(ev_halo_tree, nearest, h->halos[i].pos, range);

      sat_mass = 0;
      for (j=0; j<nearest->num_points; j++)
	if (nearest->points[j]->mvir < h->halos[i].mvir)
	  sat_mass += nearest->points[j]->mvir;
      if (sat_mass > h->halos[i].mvir) sat_mass = h->halos[i].mvir;
      if (h->halos[i].mvir)
	h->halos[i].mass_factor *= 1.0 - sat_mass / h->halos[i].mvir;
    }
    fast3tree_results_free(nearest);
  }
}

void print_evolved_halos(FILE *output, struct halo_stash *h) {
  int64_t i;
  struct tree_halo *halo;
  fprintf(output, "#ID DescID Mvir Vmax Vrms Rvir Rs Np X Y Z VX VY VZ Jx Jy Jz spin Phantom\n");
  for (i=0; i<h->num_halos; i++) {
    halo = &(h->halos[i]);

    fprintf(output, "%"PRId64" %"PRId64" %.3e %.2f %.2f %.3f %.3f %"PRId64" %.5f %.5f %.5f %.2f %.2f %.2f %.3e %.3e %.3e %.5f %"PRId64"\n",
	    halo->id, halo->descid, halo->mvir, halo->vmax, halo->vrms,
	    halo->rvir, halo->rs,
	    halo->np, halo->pos[0], halo->pos[1], halo->pos[2],
	    halo->vel[0], halo->vel[1], halo->vel[2],
	    halo->J[0], halo->J[1], halo->J[2], halo->spin,
	    halo->phantom);
  }
}

/* Assumes mvir in Msun and dt in Myr. */
inline float halo_grav_range(float mvir, float dt) {
  /* dv = dt * a = dt * (Gc*m / r^2) */
  /* So, we have dv > VEL_RESOLUTION from r=0 up to the return value: */
  return (sqrtf(fabs(mvir*Gc*dt / VEL_RESOLUTION)));
}

/* Advances simulation from a2 to a3; old a1 needed for velocity updates. */
/* special_step = 1 (first step), 0 (normal step), -1 (last step) */
void do_timestep(struct halo_stash *h, double a1, double a2, double a3, int special_step)
{
  int64_t i=0;
  float dt = (scale_to_years(a3) - scale_to_years(a2))/1.0e6; //In Myr
  float old_dt = (scale_to_years(a2) - scale_to_years(a1))/1.0e6; //In Myr
  float av_a = (a2+a3)/2.0;
  float vel_dt = dt * 1.02268944e-6 * h0 / av_a; //1 km/s to comoving Mpc/Myr/h
  float max_dist = box_size / 2.0;
  float vir_conv_factor = pow(SOFTENING_LENGTH*1.0e-3, 2);

  // 1 km/s / Myr to comoving Mpc/Myr^2/h 
  float acc_dt2 = (dt*dt / 2.0) * 1.02268944e-6 * h0 / av_a;
  float conv_const = av_a*av_a / h0 / h0;

  //Update intermediate velocities
  if (special_step != 1) { //If we're not at the first step.
#pragma omp parallel for private(i)
    for (i=0; i<h->num_halos; i++) {
      struct tree_halo *h1 = &(h->halos[i]);
      for (int64_t j=0; j<3; j++) {
	h1->vel[j] += h1->a[j]*old_dt*0.5;
	h1->a[j] = 0;
      }
    }
  }

  if (special_step<0) return; //Return if this is the final step.

  //Calculate accelerations
  fast3tree_rebuild(ev_halo_tree, h->num_halos, h->halos);
  IF_PERIODIC _fast3tree_set_minmax(ev_halo_tree, 0, box_size);
  if (special_step == 1) correct_mass_factors(h);

#pragma omp parallel default(none) shared(h,ev_halo_tree,h0,av_a,max_dist,box_size,conv_const,vir_conv_factor,dt) private(i)
  {
    struct fast3tree_results *nearest = fast3tree_results_init();
    float inv_rs, r,m, dx, dy, dz, r3, acc, rsoft2, range;

#pragma omp for
    for (i=0; i<h->num_halos; i++) {
      struct tree_halo *h1 = &(h->halos[i]);
      range = halo_grav_range(h1->mvir/h0, dt)*h0/av_a; //In comoving Mpc/h
      IF_PERIODIC {
	if (max_dist < range) range = max_dist*0.99;
	assert(fast3tree_find_sphere_periodic(ev_halo_tree, nearest, h1->pos, range));
      }
      else
	fast3tree_find_sphere(ev_halo_tree, nearest, h1->pos, range);

      inv_rs = 1000.0/h1->rs; //In comoving h/Mpc, b/c rs is in kpc/h
      for (int64_t j=0; j<nearest->num_points; j++) {
	//Calculate gravitational forces
	struct tree_halo *h2 = nearest->points[j];
#ifndef NO_PERIODIC
#define DIST(a,b) a = h1->b - h2->b; \
	if (a > max_dist) a-=box_size;		\
	else if (a < -max_dist) a+=box_size;
#else
#define DIST(a,b) a = h1->b - h2->b;
#endif
	DIST(dx,pos[0]);
	DIST(dy,pos[1]);
	DIST(dz,pos[2]);
#undef DIST
	r = sqrtf(dx*dx + dy*dy + dz*dz); // in comoving Mpc/h

	m = h1->mass_factor*ff_cached(r*inv_rs); //In Msun

	//Apply force softening:
	rsoft2 = r*r + h2->rvir*h2->rvir*vir_conv_factor;
	r3 = rsoft2 * r * conv_const; //in (real Mpc)^2 * (comoving Mpc/h)
	acc = r ? Gc*m/r3 : 0; //In km/s / Myr / (comoving Mpc/h)

#pragma omp atomic
	h2->a[0] += acc*dx;      //In km/s / Myr
#pragma omp atomic
	h2->a[1] += acc*dy;
#pragma omp atomic
	h2->a[2] += acc*dz;
      }
    }
    fast3tree_results_free(nearest);
  }

  //Update halo velocities
  if (special_step != 1) { //If not at the first step
#pragma omp parallel for private(i)
    for (i=0; i<h->num_halos; i++) {
      struct tree_halo *h1 = &(h->halos[i]);
      for (int64_t j=0; j<3; j++) h1->vel[j] += h1->a[j]*dt*0.5;
    }
  }

  //Calculate new positions
#pragma omp parallel for private(i)
  for (i=0; i<h->num_halos; i++) {
    struct tree_halo *h1 = &(h->halos[i]);
    for (int64_t j=0; j<3; j++) {
      h1->pos[j] += h1->vel[j]*vel_dt + h1->a[j]*acc_dt2;
      IF_PERIODIC {
	if (h1->pos[j] < 0) h1->pos[j] += box_size;
	if (h1->pos[j] > box_size) h1->pos[j] -= box_size;
      }
    }
  }
}
