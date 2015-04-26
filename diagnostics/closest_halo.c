#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "distance.h"
#include "kdtree.h"
#include "minimal_halo.h"
#include "masses.h"
#include "universe_time.h"
#include "gravitational_consistency.h"
#include "gravitational_consistency_vars.h"

struct min_halo *halos = NULL;
struct kdtree *halo_tree = NULL;
int num_halos = 0;
float box_size = 0; /* Automatically set; in Mpc/h */
float max_mvir = 0; /* Automatically set; in Msun */
float Om, Ol, h0; /* Omega_matter, Omega_lambda, h0 */

void find_closest_halo(void);
void load_halos(char *filename);

int main(int argc, char **argv) {
  float a1, a2, da, a;
  int i;

  if (argc < 2) {
    printf("Usage: %s halolist\n", argv[0]);
    exit(1);
  }

  /* Init cosmology and timescales */
  halo_tree = kd_create(3);
  num_halos = 0;
  max_mvir = 0;
  load_halos(argv[1]);

  find_closest_halo();
  return 0;
}

void load_halos(char *filename) {
  FILE *input;
  char buffer[1024];
  struct min_halo halo = {0};
  float max_d = 0;
  int n;
  if (!(input = fopen(filename, "r"))) {
    printf("Couldn't open file %s!\n", filename);
    exit(2);
  }
  
  halo.ax = halo.ay = halo.az = 0;
  while (fgets(buffer, 1024, input)) {
    if (buffer[0] == '#') continue;
    n = sscanf(buffer, "%d %d %f %f %f %f %d %f %f %f %f %f %f %d",&(halo.id),
	       &(halo.descendant), &(halo.mvir), &(halo.vmax), &(halo.rvir),
	       &(halo.rs), &(halo.np), &(halo.x), &(halo.y), &(halo.z),
	       &(halo.vx), &(halo.vy), &(halo.vz), &(halo.phantom_id));
    if (n < 13) continue;
    if (n < 14) halo.phantom_id = 0;
    if (!(halo.mvir > 0)) continue;
    if (!(halo.rvir > 0)) continue;
    if (!(halo.vmax > 0)) continue;
    halo.vmax = log10(halo.vmax);

    if (!(num_halos % 1000)) {
      halos = (struct min_halo *)
	realloc(halos, sizeof(struct min_halo)*(num_halos+1000));
      if (!halos) {
	printf("Out of memory trying to load halos!\n");
	exit(1);
      }
    }
    halos[num_halos] = halo;
    num_halos++;

    if (max_mvir < halo.mvir) max_mvir = halo.mvir;
    if (max_d < halo.x) max_d = halo.x;
  }
  fclose(input);
  box_size = (int)(max_d + 0.5); //In Mpc/h comoving coordinates.
  max_mvir /= h0; //Now in Msun
}


/* Max range is the maximum distance over which the gravitational effects of
   the largest halo will matter. */
void build_halo_tree(float max_range) {
  int i,j,k,l,m,n;
  float nx, ny, nz;
  if (max_range > box_size/4) { max_range = box_size/4; }

  float inv_max_range = box_size - max_range;

  if (inv_max_range < max_range) { //Weird things will happen, so exit
    printf("Velocity resolution is too precise (would need to consider too\n"
	   "many halos).  Reduce velocity resolution and recompile.\n");
    exit(2);
  }

  /* Since the kdtree library doesn't have cyclic boundary conditions,
     we need to duplicate a few halos across the boundaries of the box. */
  for (i=0; i<num_halos; i++) {
    j = (halos[i].x<max_range || halos[i].x>inv_max_range) ? 1 : 0;
    k = (halos[i].y<max_range || halos[i].y>inv_max_range) ? 1 : 0;
    l = (halos[i].z<max_range || halos[i].z>inv_max_range) ? 1 : 0;
    if (j||k||l) {
      nx = halos[i].x + box_size * ((halos[i].x<max_range) ? 1 : -1);
      ny = halos[i].y + box_size * ((halos[i].y<max_range) ? 1 : -1);
      nz = halos[i].z + box_size * ((halos[i].z<max_range) ? 1 : -1);
      for (; j>=0; j--)
	for (m=k; m>=0; m--)
	  for (n=l; n>=0; n--)
	    kd_insert3f(halo_tree, (j ? nx : halos[i].x), 
			(m ? ny : halos[i].y), (n ? nz : halos[i].z), 
			&(halos[i]));
    }
    else {
      kd_insert3f(halo_tree, halos[i].x, halos[i].y, halos[i].z, &(halos[i]));
    }
  }
}

void destroy_halo_tree(void) {
  kd_clear(halo_tree);
}


/* Advances simulation from a2 to a3; old a1 needed for velocity updates. */
/* special_step = 1 (first step), 0 (normal step), -1 (last step) */
void find_closest_halo(void)
{
  int i,j;
  struct kdres *nearest;
  struct min_halo *h1, *h2, *md;
  float min_d = box_size;
  float max_dist = box_size/2.0;
  float range, dx, dy, dz, r;

  build_halo_tree(20);
  for (i=0; i<num_halos; i++) {
    h1 = &(halos[i]);
    range = h1->rvir / 1000.0;
    md = 0;
    min_d = box_size;
    while ((!md) && (range < (h1->rvir*30.0/1000.0))) {
      nearest = kd_nearest_range3f(halo_tree, h1->x, h1->y, h1->z, range);
      if (!kd_res_size(nearest)) {
	kd_res_free(nearest);
	range *= 2;
	continue; //Skip if no close halos.
      }

      do {
	//Calculate gravitational forces
	h2 = (struct min_halo *)kd_res_item_data(nearest);
#define DIST(a,b) a = h1->b - h2->b; \
	if (a > max_dist) a-=box_size;		\
	else if (a < -max_dist) a+=box_size;
	DIST(dx,x);
	DIST(dy,y);
	DIST(dz,z);
#undef DIST
	r = sqrtf(dx*dx + dy*dy + dz*dz); // in comoving Mpc/h
	if ((h2->id != h1->id) && (r < min_d) && 
	    (fabs(h2->vmax - h1->vmax) < VMAX_LOG10_TOL)) {
	  min_d = r;
	  md = h2;
	}
      } while (kd_res_next(nearest));
      kd_res_free(nearest);
      range *=2;
    }
    printf("%d %e %f %f ", h1->id, h1->mvir, h1->rvir, pow(10, h1->vmax));
    if (md) {
      printf("%d %f\n",  md->id, min_d);
    } else {
      printf("%d %f\n",  -1, range/2.0);
    }
  }
  destroy_halo_tree();
}
