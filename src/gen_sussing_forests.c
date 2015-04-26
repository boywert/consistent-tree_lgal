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

#define forest_id tidal_id
#define CRITICAL_DENSITY 2.77459457e11 // 3H^2/8piG in (Msun / h) / (Mpc / h)^3

//#define FIX_BOLSHOI_SPIN


struct halo_stash now={0}, evolved = {0};
int64_t *forest_list = NULL;
int64_t num_forests = 0;
float box_size=0;
float max_mvir=0;
float min_mvir=0;
int64_t children = 0;

void write_ahf_halos(struct halo_stash *h, int64_t onum, float scale);
void wait_for_children(int all);
int sort_by_forest(const void *a, const void *b);
void load_forests(struct halo_stash *h);

struct id_num {
  int64_t id, num;
};

int main(int argc, char **argv) {
  int64_t i, j, num_outputs=0;
  float *output_scales=NULL;
  int64_t *outputs=NULL;
  char buffer[1024];
  int64_t full_halo_count=0;

  if (argc==1) {
    fprintf(stderr, "Consistent Trees, Version %s\n", TREE_VERSION);
    fprintf(stderr, "%s.  See the LICENSE file for redistribution details.\n", TREE_COPYRIGHT);
    fprintf(stderr, "Usage: %s options.cfg\n", argv[0]); exit(1);
  }
  if (argc>1) grav_config(argv[1], 1);
  read_outputs(&output_scales, &outputs, &num_outputs);
  gen_ff_cache();
  clear_halo_stash(&now);

  snprintf(buffer, 1024, "%s/really_consistent_%"PRId64".list", OUTBASE, outputs[num_outputs-1]);
  load_halos(buffer, &now, output_scales[num_outputs-1], 0);
  build_id_conv_list(&now);
  load_forests(&now);
  qsort(now.halos, now.num_halos, sizeof(struct tree_halo), sort_by_forest);
  build_id_conv_list(&now);
  write_ahf_halos(&now, num_outputs-1, output_scales[num_outputs-1]);

  for (i = num_outputs-1; i>=0; i--) {
    clear_halo_stash(&evolved);
    evolved = now;
    zero_halo_stash(&now);
    
    if (i>0) {
      snprintf(buffer, 1024, "%s/really_consistent_%"PRId64".list", OUTBASE, outputs[i-1]);
      load_halos(buffer, &now, output_scales[i-1], 0);
      build_id_conv_list(&now);
      for (j=0; j<now.num_halos; j++) {
	int64_t desc_index = id_to_index(evolved, now.halos[j].descid);
	assert(desc_index > -1);
	now.halos[j].forest_id = evolved.halos[desc_index].forest_id;
      }
    
      qsort(now.halos, now.num_halos, sizeof(struct tree_halo), sort_by_forest);
      build_id_conv_list(&now);
      write_ahf_halos(&now, i-1, output_scales[i-1]);
    }

    snprintf(buffer, 1024, "%s/sussing_%"PRId64".list", TREE_OUTBASE, i);
    FILE *output = check_fopen(buffer, "w");

    int64_t k=0, last_k=0, l=0, last_l=0;
    for (j=0; j<num_forests; j++) {
      struct id_num forest;
      forest.id = forest_list[j];
      if (last_k < evolved.num_halos)
	for (; last_k<evolved.num_halos; last_k++)
	  if (evolved.halos[last_k].forest_id != forest.id) break;
      forest.num = last_k - k;
      fwrite(&forest, sizeof(struct id_num), 1, output);
      if (!forest.num) continue;

      for (; k<last_k; k++) {
	struct id_num halo;
	halo.id = evolved.halos[k].id;
	if (last_l < now.num_halos)
	  for (; last_l < now.num_halos; last_l++)
	    if (now.halos[last_l].descid != halo.id) break;
	halo.num = last_l - l;
	fwrite(&halo, sizeof(struct id_num), 1, output);
	if (!halo.num) continue;
	for (; l<last_l; l++)
	  fwrite(&(now.halos[l].id), sizeof(int64_t), 1, output);
      }
    }
    fclose(output);
  }

  FILE **inputs = check_realloc(NULL, num_outputs*sizeof(FILE *), "Fh");
  for (i=0; i<num_outputs; i++) {
    snprintf(buffer, 1024, "%s/sussing_%"PRId64".list", TREE_OUTBASE, i);
    inputs[i] = check_fopen(buffer, "r");
  }

  snprintf(buffer, 1024, "%s/sussing_tree.list", TREE_OUTBASE);
  FILE *tree = check_fopen(buffer, "w");
  fprintf(tree, "1\nConsistent Trees version %s\n", TREE_VERSION);  
  int64_t size_offset = ftello(tree);
  fprintf(tree, "XXXXXXXXXXXXXXXXXX\n");

  snprintf(buffer, 1024, "%s/sussing_forests.list", TREE_OUTBASE);
  FILE *fcounts = check_fopen(buffer, "w");
  fprintf(fcounts, "#ForestID Num_Halos_Total Halos_In_Snap_0 Halos_In_Snap_1 ...\n");

  int64_t *halo_counts = check_realloc(NULL, num_outputs*sizeof(int64_t), "Halo Counts");
  for (j=0; j<num_forests; j++) {
    int64_t fid = forest_list[j];
    int64_t total_halos = 0;
    memset(halo_counts, 0, sizeof(int64_t)*num_outputs);
    for (i=0; i<num_outputs; i++) {
      struct id_num forest;
      fread(&forest, sizeof(struct id_num), 1, inputs[i]);
      assert(forest.id == fid);
      halo_counts[i] = forest.num;
      total_halos += forest.num;
      full_halo_count += forest.num;
      int64_t k, l;
      for (k=0; k<forest.num; k++) {
	struct id_num halo;
	fread(&halo, sizeof(struct id_num), 1, inputs[i]);
	fprintf(tree, "%"PRId64" %"PRId64"\n", halo.id, halo.num);
	for (l=0; l<halo.num; l++) {
	  int64_t progid;
	  fread(&progid, sizeof(int64_t), 1, inputs[i]);
	  fprintf(tree, "%"PRId64"\n", progid);
	}
      }
    }

    fprintf(fcounts, "%"PRId64" %"PRId64, fid, total_halos);
    for (i=0; i<num_outputs; i++) fprintf(fcounts, " %"PRId64, halo_counts[i]);
    fprintf(fcounts, "\n");
  }
  fprintf(tree, "END\n");
  fseeko(tree, size_offset, SEEK_SET);
  fprintf(tree, "%-18"PRId64"\n", full_halo_count);
  fclose(tree);
  fclose(fcounts);

  for (i=0; i<num_outputs; i++) {
    fclose(inputs[i]);
    snprintf(buffer, 1024, "%s/sussing_%"PRId64".list", TREE_OUTBASE, i);
    unlink(buffer);
  }
  return 0;
}


float dens_200c(float z) {
  float inv_vol = pow(1.0+z, 3);
  float dens_200m = 200.0*(Om*inv_vol+Ol)/(Om*inv_vol);
  float background = CRITICAL_DENSITY*Om;
  return dens_200m*background;
}
 
void write_ahf_halos(struct halo_stash *h, int64_t onum, float scale) {
  int64_t *num_subs = check_realloc(NULL, sizeof(int64_t)*h->num_halos, "Subcounts");
  memset(num_subs, 0, sizeof(int64_t)*h->num_halos);
  int64_t n,i,j;
  for (i=0; i<h->num_halos; i++) {
    int64_t upid = i;
    while (h->halos[upid].upid > -1) upid = id_to_index(*h, h->halos[upid].upid);
    num_subs[upid]++;
  }

  char buffer[1024];
  snprintf(buffer, 1024, "%s/sussing_%03"PRId64".z%.3f.AHF_halos", TREE_OUTBASE, onum, 1.0/scale-1.0);
  FILE *output = check_fopen(buffer, "w");
  fprintf(output, "#ID(1)  hostHalo(2)     numSubStruct(3) M200c(4) npart(5)        Xc(6)   Yc(7)   Zc(8)   VXc(9)  VYc(10) VZc(11) R200c(12)        Rmax(13)        r2(14)  mbp_offset(15)  com_offset(16)  Vmax(17)        v_esc(18)       sigV(19)        lambda(20)      lambdaE(21)     Lx(22)  Ly(23)  Lz(24) cNFW(25)\n");

  float dens = dens_200c(1.0/scale - 1.0)*(4.0*M_PI/3.0);
  //SUSSING_MASS_FIELD = 3;
  for (n=0; n<h->num_halos; n++) {
    i = n; //sort_order[n];
    struct tree_halo *th = h->halos+i;
    int64_t upid = i;
    while (h->halos[upid].upid > -1) upid = id_to_index(*h, h->halos[upid].upid);
    upid = (upid == i) ? 0 : h->halos[upid].id;
    float L[3] = {0}, Ln=0;
    float mass = th->mvir;
    float radius = th->rvir;
    if (SUSSING_MASS_FIELD > -1) {
      mass = th->extra_params[SUSSING_MASS_FIELD];
      if (!(mass>0)) {
	convert_mvir_to_delta(th->mvir, th->rvir, th->rs, scale, 200, 'c', &mass, &radius);
      }
      radius = cbrt(mass/dens)*1e3;
      
    }
    float spin = th->spin;
    for (j=0; j<3; j++) Ln+=th->J[j]*th->J[j];
    Ln = sqrt(Ln);
    if (Ln) for (j=0; j<3; j++) L[j] = th->J[j]/Ln;

#ifdef FIX_BOLSHOI_SPIN
    float t_u = th->extra_params[19];
    if (fabs(1.0-t_u) != 0)
      spin *= sqrt(fabs(1.0-t_u*2.0)) / sqrt(fabs(1.0-t_u));
#endif
    float cNFW = 10;
    if (th->rs > 0) cNFW = radius / th->rs;
    fprintf(output, "%"PRId64" %"PRId64" %"PRId64" %e %"PRId64" %.6f %.6f %.6f %f %f %f %f %f 2e38 2e38 2e38 %f 2e38 %f 2e38 %f %f %f %f %f\n", th->id, upid, num_subs[i], mass, th->np, th->pos[0]*1e3, th->pos[1]*1e3, th->pos[2]*1e3, th->vel[0], th->vel[1], th->vel[2], radius, th->rs*2.1626, th->vmax, th->vrms, spin, L[0], L[1], L[2], cNFW);
  }
  fclose(output);
  free(num_subs);
}

int sort_ascending(const void *a, const void *b) {
  const int64_t *c = a;
  const int64_t *d = b;
  if (*c < *d) return -1;
  return 1;
}

void load_forests(struct halo_stash *h) {
  int64_t key, val;
  char buffer[1024];
  snprintf(buffer, 1024, "%s/forests.list", TREE_OUTBASE);
  FILE *input = check_fopen(buffer, "r");
  fgets(buffer, 1024, input); //Header
  struct inthash *forests = new_inthash();
  for (key=0; key<h->num_halos; key++)
    h->halos[key].forest_id = -1;

  while (fgets(buffer, 1024, input)) {
    if (sscanf(buffer, "%"SCNd64" %"SCNd64, &key, &val) != 2) continue;
    int64_t idx = id_to_index(*h, key);
    assert(idx > -1);
    h->halos[idx].forest_id = val;
    ih_setint64(forests, val, 1);
  }
  fclose(input);

  for (key=0; key<h->num_halos; key++)
    assert(h->halos[key].forest_id != -1);

  forest_list = ih_keylist(forests);
  num_forests = forests->elems;
  free_inthash(forests);
  qsort(forest_list, num_forests, sizeof(int64_t), sort_ascending);
}

inline int64_t id_to_index(struct halo_stash h, int64_t id) {
  if (id < h.min_id || id > h.max_id) return -1;
  return (h.id_conv[id-h.min_id]);
}

int sort_by_forest(const void *a, const void *b) {
  const struct tree_halo *c = a; //now.halos + *((int64_t *)a);
  const struct tree_halo *d = b; //now.halos + *((int64_t *)b);
  if (c->forest_id < d->forest_id) return -1;
  if (c->forest_id > d->forest_id) return 1;
  if (c->descid > -1) {
    assert(d->descid > -1);
    int64_t c_idx = id_to_index(evolved, c->descid);
    int64_t d_idx = id_to_index(evolved, d->descid);
    if (c_idx < d_idx) return -1;
    if (c_idx > d_idx) return 1;
  }
  if (c->mvir > d->mvir) return -1;
  if (c->mvir < d->mvir) return 1;
  if (c->vmax > d->vmax) return -1;
  if (c->vmax < d->vmax) return 1;
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

