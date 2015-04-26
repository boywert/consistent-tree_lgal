#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include "gravitational_consistency.h"
#include "gravitational_consistency_vars.h"
#include "masses.h"
#include "distance.h"
#include "universe_time.h"
#include "tidal_lib.h"
#include "stringparse.h"
#include "halo_io.h"
#include "grav_config.h"
#include "check_syscalls.h"

struct halo_stash now = {0};
struct halo_stash evolved = {0};
float box_size = 0; /* Automatically set; in Mpc/h */
float max_mvir = 0; /* Automatically set; in Msun */
float min_mvir = 0;

void build_mmp_list(void);
void print_tidal_forces(int64_t output_num, float a, struct halo_stash *h, struct halo_stash *evolved);

int main(int argc, char **argv) {
  int64_t i, num;
  int64_t num_outputs=0;
  float *output_scales=NULL;
  int64_t *outputs=NULL;
  char buffer[1024];

  if (argc < 3) {
    printf("Usage: %s configfile output_num\n", argv[0]);
    exit(1);
  }

  grav_config(argv[1], 0);
  num = atoi(argv[2]);
  TIDAL_FORCE_LIMIT = 0.0001;

  init_cosmology(Om, Ol, h0);
  init_time_table(Om, h0);
  read_outputs(&output_scales, &outputs, &num_outputs);
  gen_ff_cache();

  for (i=0; i<num_outputs; i++)
    if (outputs[i]==num) break;

  if (i==num_outputs) {
    printf("Couldn't find output %"PRId64"!\n", num);
    exit(1);
  }
  if (i==num_outputs-1) {
    printf("Not meaningful to calculate tidal fields for last output!\n");
    exit(1);
  }

  max_mvir = 0;
  snprintf(buffer, 1024, "%s/out_%"PRId64".list", INBASE, outputs[i]);
  load_halos(buffer, &now, output_scales[i], 0);
  build_id_conv_list(&now);
  snprintf(buffer, 1024, "%s/out_%"PRId64".list", INBASE, outputs[i+1]);
  load_halos(buffer, &evolved, output_scales[i+1], 0);
  build_id_conv_list(&evolved);
  build_mmp_list();
  calc_tidal_forces(&now, output_scales[i], output_scales[i+1]);
  print_tidal_forces(outputs[i], output_scales[i], &now, &evolved);
  
  return 0;
}


inline int64_t id_to_index(struct halo_stash h, int64_t id) {
  if (id < h.min_id || id > h.max_id) return -1;
  return (h.id_conv[id-h.min_id]);
}

void build_mmp_list(void)
{
  int64_t j, eindex, index;
  //Build up the most-massive progenitor list
  //This simultaneously tags halos w/o progenitors
  for (j=0; j<evolved.num_halos; j++) evolved.halos[j].mmp_id = -1;

  for (j=0; j<now.num_halos; j++) {
    now.halos[j].flags -= (now.halos[j].flags & MMP_FLAG);
    if (now.halos[j].flags & DEAD_HALO_FLAG) {
      now.halos[j].flags |= MMP_FLAG;
      continue;
    }
 
    eindex = id_to_index(evolved, now.halos[j].descid);
    if (eindex < 0) { //Descendant halo no longer exists
      now.halos[j].flags |= MMP_FLAG | DEAD_HALO_FLAG;
      continue;
    }

    index = id_to_index(now, evolved.halos[eindex].mmp_id);
    if (index < 0) { //If no MMP already identified.
      evolved.halos[eindex].mmp_id = now.halos[j].id;
      now.halos[j].flags |= MMP_FLAG;
      continue;
    }

    //Otherwise, check to see if we're more massive than the previous MMP.
    if (now.halos[j].mvir >= now.halos[index].mvir) {
      evolved.halos[eindex].mmp_id = now.halos[j].id;
      now.halos[index].flags -= (now.halos[index].flags & MMP_FLAG); //Safer than &= ~MMP_FLAG
      now.halos[j].flags |= MMP_FLAG;
    }
  }
}


void print_tidal_forces(int64_t output_num, float a, struct halo_stash *h, struct halo_stash *evolved) {
  char buffer[1024];
  FILE *mmp_out, *nmmp_out, *output;
  int64_t mmp, j, eindex, tidal_desc_id, nindex, bfnmmp;
  float tidal_mvir;
  snprintf(buffer, 1024, "%s/extended_tidal_forces_mmp_%"PRId64".list", OUTBASE, output_num);
  mmp_out = check_fopen(buffer, "w");
  fprintf(mmp_out, "#ID DescID TidalId TidalDescId TidalForce Rvir Mvir TidalMvir Vmax\n");
  fprintf(mmp_out, "#Scale factor: %f\n", a);
  snprintf(buffer, 1024, "%s/extended_tidal_forces_nmmp_%"PRId64".list", OUTBASE, output_num);
  nmmp_out = check_fopen(buffer, "w");
  fprintf(nmmp_out, "#ID DescID TidalId TidalDescId TidalForce Rvir Mvir TidalMvir Vmax Confirmed_Nmmp\n");
  fprintf(nmmp_out, "#Scale factor: %f\n", a);
  for (j=0; j<h->num_halos; j++) {
    if (h->halos[j].flags & DEAD_HALO_FLAG)
      continue;
 
    eindex = id_to_index(*evolved, h->halos[j].descid);
    if (eindex < 0) //Descendant halo no longer exists
      continue;
    bfnmmp = 1;
    mmp = (evolved->halos[eindex].mmp_id == h->halos[j].id) ? 1 : 0;
    if (mmp && really_beyond_mmp_ratio(h->halos[j], evolved->halos[eindex])) {
      mmp = 0;
      bfnmmp = 0;
    }
    output = (mmp) ? mmp_out : nmmp_out;

    nindex = id_to_index(*h, h->halos[j].tidal_id);
    if (nindex>-1) {
      tidal_desc_id = h->halos[nindex].descid;
      tidal_mvir = h->halos[nindex].mvir;
    }
    else { tidal_desc_id = -1; tidal_mvir = 0; }

    fprintf(output, "%"PRId64" %"PRId64" %"PRId64" %"PRId64" %.3e %f %.3e %.3e %f %"PRId64"\n",
            h->halos[j].id, h->halos[j].descid, h->halos[j].tidal_id,
            tidal_desc_id, h->halos[j].tidal_force, h->halos[j].rvir, 
            h->halos[j].mvir, tidal_mvir, h->halos[j].vmax, bfnmmp);
  }
  fclose(mmp_out);
  fclose(nmmp_out);
}
