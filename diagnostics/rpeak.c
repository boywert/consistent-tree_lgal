#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <sys/stat.h>

#define EXTRA_HALO_INFO struct halo *vpeak_halo, *mpeak_halo, *acc_halo; int64_t sub_ok; int64_t flyby; struct halo *mhalo[5]; int64_t incomplete; int64_t unresolved_mergers;
#include "../read_tree/read_tree.h"
#include "../read_tree/read_tree.c"
#include "../read_tree/check_syscalls.h"
#include "../src/universe_time.h"
#include "rpeak.h"

float box_size = 250;
void process_tree(void);
float calc_radius(struct halo *h, struct halo *p);
void calc_host_mvir_rvir(struct halo *h, struct halo *p, float *mvir_acc, float *rvir_acc);
struct halo *find_parent_at_scale(float scale, struct halo *h);

#define MVIR_CUT 2e10

struct halo_data *hd = NULL;
int64_t num_hd = 0;

struct thalo {
  float pos[3], r, m;
};
#define FAST3TREE_TYPE struct thalo 
#include "../src/fast3tree.c"

struct thalo *th = NULL;
int64_t *num_th = NULL;
int64_t total_th = 0;
struct fast3tree **trees = NULL;
int64_t num_trees = 0;
double *tree_scales = NULL;
struct fast3tree_results *res = NULL;
float RATIOS[5] = {40,20,10,5,2.5};

void load_kdtrees(char *filename);
void find_nearest_larger_halo(struct halo *h, float *nr, float *nm);

int main(int argc, char **argv) {
  int64_t i;
  FILE *output;
  if (argc<4) {
    printf("Usage: %s box_size AllLC.dat tree.dat ...\n", argv[0]);
    exit(1);
  }

  box_size = atof(argv[1]);
  init_time_table(0.27, 0.7);
  //  load_kdtrees(argv[2]);
  
  printf("#Scale Id Mvir_acc(Host) Rvir_acc(Host) Mnow Vnow Mpeak Vpeak RVpeak RMpeak RVpeak/Rvir Vacc Mistracked? VPflyby VPissub MPflyby MPissub VPdif_Sub MPdif_Sub Vpeak_Scale D_to_Larger_at_Vpeak_over_Rvir M_Larger_at_Vpeak R_Last1:%g R_Last1:%g R_Last1:%g R_Last1:%g R_Last1:%g Unresolved_mergers\n", RATIOS[0], RATIOS[1], RATIOS[2],
RATIOS[3], RATIOS[4]);
  output = fopen("host_flyby.dat", "w");
  fclose(output);

  for (i=3; i<argc; i++) {
    read_tree(argv[i]);
    process_tree();
    delete_tree();
  }
  return 0;
}

void calc_last_mergers(void) {
  int64_t i, j;
  for (i=0; i<all_halos.num_halos; i++) {
    struct halo *h = &(all_halos.halos[i]);
    for (j=0; j<5; j++) h->mhalo[j] = NULL;
    h->unresolved_mergers = 0;
  }

  for (i=0; i<all_halos.num_halos; i++) {
    struct halo *h = &(all_halos.halos[i]);
    if (h->pid < 0 || !h->parent || !h->prog) continue;
    if (h->prog->parent) continue;
    for (j=0; j<5; j++) {
      if (!h->parent->mhalo[j] && h->prog->mvir*RATIOS[j] > h->parent->mvir) {
	struct halo *desc = h->parent;
	while (desc && (!desc->mhalo[j] || 
			desc->mhalo[j]->scale < h->parent->scale)) {
	  desc->mhalo[j] = h->parent;
	  if (!desc->mmp) break;
	  desc = desc->desc;
	}
      }
    }
  }

  for (i=0; i<all_halos.num_halos; i++) {
    struct halo *h = &(all_halos.halos[i]);
    for (j=0; j<5; j++) {
      struct halo *prog = h;
      if (h->mhalo[j]) continue;
      h->unresolved_mergers = 1;
      while (prog->prog) prog = prog->prog;
      h->mhalo[j] = prog;
    }
  }
}

void calc_mass_vmax_peak(struct halo *h) {
  if (h->mpeak_halo) return;
  if (h->prog) {
    calc_mass_vmax_peak(h->prog);
    h->flyby = h->prog->flyby;
    if (h->vmax > h->prog->vpeak_halo->vmax) h->vpeak_halo = h;
    else h->vpeak_halo = h->prog->vpeak_halo;

    if (h->mvir > h->prog->mpeak_halo->mvir) h->mpeak_halo = h;
    else h->mpeak_halo = h->prog->mpeak_halo;

    if (h->upid>-1) { // If we are a subhalo
      h->acc_halo = h->prog->acc_halo;
      h->flyby = 1;
    } else {
      h->acc_halo = h;
    }
  } else {
    h->acc_halo = h->mpeak_halo = h->vpeak_halo = h;    
    h->flyby = (h->pid < 0) ? 0 : 1;
  }

  if (h->acc_halo->pid < 0) h->sub_ok = 1;
  else h->sub_ok = 0;
}

void check_if_near_boundary(void) {
  int64_t i, max_i=0, j;
  float max_scale = 0;
  for (i=0; i<halo_tree.num_lists; i++) {
    if (halo_tree.halo_lists[i].scale > max_scale) {
      max_scale = halo_tree.halo_lists[i].scale;
      max_i = i;
    }
  }

  for (j=0; j<3; j++) {
    float min = 1e10;
    float max = 0;
    for (i=0; i<halo_tree.halo_lists[max_i].num_halos; i++) {
      struct halo *h = halo_tree.halo_lists[max_i].halos+i;
      float pos = h->pos[j];
      if (!j) h->incomplete = 0;
      if (pos < min) min = pos;
      if (pos > max) max = pos;
    }
    for (i=0; i<halo_tree.halo_lists[max_i].num_halos; i++) {
      struct halo *h = halo_tree.halo_lists[max_i].halos+i;
      float pos = h->pos[j];
      if ((pos < min+2.0) || (pos > max-2.0)) {
	while (h) { h->incomplete = 1; h = h->prog; }
      }
    }
  }
}


void process_tree(void) {
  float rvpeak, rmpeak, mvir_acc, rvir_acc, ratio;
  int64_t i, j;
  struct halo *h;
  FILE *flyby;
  flyby = fopen("host_flyby.dat", "a");
  calc_last_mergers();
  check_if_near_boundary();
  for (i=0; i<all_halos.num_halos; i++) {
    h = &(all_halos.halos[i]);
    calc_mass_vmax_peak(h);
    if (h->pid < 0) {
      fprintf(flyby, "%.5f %"PRId64" %.3e %.4f %"PRId64"\n", h->scale, h->id, h->mpeak_halo->mvir, h->vpeak_halo->vmax, h->flyby);
      continue;
    }
    if (!h->parent) continue;
    if (h->incomplete) continue;
    
    rvpeak = calc_radius(h->vpeak_halo, h->parent);
    rmpeak = calc_radius(h->mpeak_halo, h->parent);
    calc_host_mvir_rvir(h->acc_halo, h->parent, &mvir_acc, &rvir_acc);
    //calc_host_mvir_rvir(h->vpeak_halo, h->parent, &mvir_acc, &rvir_acc);
    int64_t vpeak_issub = (h->vpeak_halo->pid < 0) ? 0 : 1;
    int64_t mpeak_issub = (h->mpeak_halo->pid < 0) ? 0 : 1;

    int64_t mp_difsub=0, vp_difsub=0;
    struct halo *p = find_parent_at_scale(h->vpeak_halo->scale, h->parent);
    if (vpeak_issub && p && (p->id != h->vpeak_halo->pid) && (p->id != h->vpeak_halo->upid)) vp_difsub = 1;
    p = find_parent_at_scale(h->mpeak_halo->scale, h->parent);  
    if (mpeak_issub && p && (p->id != h->mpeak_halo->pid) && (p->id != h->mpeak_halo->upid)) mp_difsub = 1;
    float d_lh=0, m_lh=0;
    float rm[5];
    for (j=0; j<5; j++) rm[j] = calc_radius(h->mhalo[j], h->parent);

    //find_nearest_larger_halo(h->vpeak_halo, &d_lh, &m_lh);

    ratio = (rvir_acc > 0) ? (rvpeak / rvir_acc) : 0.0;
    printf("%.5f %"PRId64" %.3e %.4f %.3e %.2f %.3e %.2f %f %f %f %f %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %f %f %e %f %f %f %f %f %"PRId64"\n", h->scale, h->id, mvir_acc, rvir_acc, h->mvir, h->vmax, h->mpeak_halo->mvir, h->vpeak_halo->vmax, rvpeak, rmpeak, ratio, h->acc_halo->vmax, (1-h->sub_ok), h->vpeak_halo->flyby, vpeak_issub, h->mpeak_halo->flyby, mpeak_issub, vp_difsub, mp_difsub, h->vpeak_halo->scale, d_lh, m_lh, rm[0], rm[1], rm[2], rm[3], rm[4], h->unresolved_mergers);
  }
  fclose(flyby);
}



float calc_dist(struct halo *h1, struct halo *h2) {
  float sum = 0, dt, av_a, cmvng_mpc_per_kms, dx;
  int64_t i;
  struct halo h3;
  memcpy(&h3, h2, sizeof(struct halo));
  h2 = &h3;
  if (h2->scale != h1->scale) {
    dt = scale_to_years(h2->scale) - scale_to_years(h1->scale);
    av_a = (h1->scale+h2->scale)/2.0;
    cmvng_mpc_per_kms = 1.022e-6/av_a*dt/1.0e6;
    for (i=0; i<3; i++)
      h2->pos[i] -= cmvng_mpc_per_kms*h2->vel[i];
  }

  for (i=0; i<3; i++) {
    dx = fabs(h1->pos[i] - h2->pos[i]);
    if (dx > box_size / 2.0) dx = box_size - dx;
    sum += dx * dx;
  }
  return sqrt(sum);
}

struct halo* find_earliest_parent(struct halo *h) {
  struct halo *p = h;
  while (p->prog) p = p->prog;
  return p;
}

struct halo *find_parent_at_scale(float scale, struct halo *h) {
  struct halo *p = h;
  while (p->scale > scale) {
    p = p->prog;
    if (!p) return 0;
  }
  return p;
}

float calc_radius(struct halo *h, struct halo *p) {
  struct halo *p2 = find_parent_at_scale(h->scale, p);
  if (!p2) p2 = find_earliest_parent(p);
  if (!p2) return 0;
  return calc_dist(h, p2);
}

void calc_host_mvir_rvir(struct halo *h, struct halo *p, float *mvir_acc, float *rvir_acc) {
 struct halo *p2 = find_parent_at_scale(h->scale, p);
  if (!p2) p2 = find_earliest_parent(p);
  if (!p2) { *mvir_acc = *rvir_acc = 0; return; }
  *mvir_acc = p2->mvir;
  *rvir_acc = p2->rvir/1.0e3;
}

void load_bin_cat(char *filename, double scale) {
  FILE *input = check_fopen(filename, "r");
  struct stat file_stat;
  int64_t cat_halos = 0;
  fstat(fileno(input), &file_stat);
  if (file_stat.st_size % sizeof(struct halo_data)) {
    fprintf(stderr, "Size of file %s is not divisible by structure size (%"PRId64")!\n", filename, (int64_t)(sizeof(struct halo_data)));
    exit(1);
  }
  
  cat_halos = file_stat.st_size /  sizeof(struct halo_data);
  hd = (struct halo_data *)
    check_realloc(hd, sizeof(struct halo_data)*(cat_halos),
		  "Allocating halo data");

  //Read everything into the halos array.
  fread(hd, sizeof(struct halo_data), cat_halos, input);
  fclose(input);

  th = check_realloc(th, sizeof(struct thalo)*(total_th+cat_halos),
		     "Allocating tree halo data");
  tree_scales = check_realloc(tree_scales, sizeof(double)*(num_trees+1),
		     "Allocating tree sizes");
  num_th = check_realloc(num_th, sizeof(int64_t)*(num_trees+1),
		     "Allocating tree sizes");
  num_th[num_trees] = 0;
  tree_scales[num_trees] = scale;
  for (int64_t i=0; i<cat_halos; i++) {
    if (hd[i].mass > MVIR_CUT) {
      memcpy(th[total_th].pos, hd[i].p, sizeof(float)*3);
      th[total_th].r = hd[i].rvir;
      th[total_th].m = hd[i].mass;
      total_th++;
      num_th[num_trees]++;
    }
  }
  num_trees++;
}


void create_kdtrees(void) {
  int64_t i, t_start=0;
  trees = check_realloc(NULL, sizeof(struct fast3tree *)*num_trees,
			"Allocating trees");
  for (i=0; i<num_trees; i++) {
    trees[i] = fast3tree_init(num_th[i], th+t_start);
    _fast3tree_set_minmax(trees[i], 0, box_size);
    t_start += num_th[i];
  }
  res = fast3tree_results_init();
}

float dist(float *p1, float *p2) {
  int64_t j;
  double ds=0, dx;
  for (j=0; j<3; j++) {
    dx = fabs(p1[j]-p2[j]);
    if (dx > box_size/2.0) dx = box_size-dx;
    ds+=dx*dx;
  }
  return sqrt(ds);
}

void find_nearest_larger_halo(struct halo *h, float *nr, float *nm) {
  *nr = *nm = 0;
  if (h->scale>1.0 || h->scale < 0) return;
  int64_t i=0;
  for (i=0; i<num_trees-1; i++)
    if (fabs(tree_scales[i+1]-h->scale) > fabs(tree_scales[i]-h->scale)) break;

  double d2=1e4;
  float r = h->rvir;
  for (; r<1e3*box_size/2.0; r*=2.0) {
    fast3tree_find_sphere_periodic(trees[i], res, h->pos, r/1.0e3);
    for (int64_t j=0; j<res->num_points; j++) {
      struct thalo *dh = res->points[j];
      if (dh->r <= h->rvir*1.02) continue;
      double d = dist(h->pos, dh->pos)*1.0e3;
      if (d/dh->r < d2) {
	d2 = d/dh->r;
	*nm = dh->m;
	*nr = d2;
      }
    }
    if (r/3.0e3 > d2) break;
  }
}

float readfloat(FILE *input) {
  char buffer[50] = {0};
  float result;
  fgets(buffer, 50, input);
  sscanf(buffer, "%f", &result);
  return result;
}

void load_kdtrees(char *filename)
{
  float scale;
  FILE *input = check_fopen(filename, "r");
  char buffer[1024], halo_file[1024];
  float Om, Ol, h0; 

  Om = readfloat(input);
  Ol = readfloat(input);
  h0 = readfloat(input);
  if (Om <= 0 || Om > 1 || Ol <= 0 || Ol > 1 ||
      h0 <= 0 || h0 > 1) {
    printf("Incorrect format for halo catalog in %s!\n", filename);
    exit(1);
  }
  if (fabs(Om+Ol-1) > 0.01) {
    printf("Flat cosmologies only.\n");
    exit(4);
  }

  init_time_table(Om, h0);
  while (fgets(buffer, 1024, input)) {
    if (!strlen(buffer)) continue;
    if (sscanf(buffer, "%f %[^\n]", &scale, halo_file) != 2) continue;
    load_bin_cat(halo_file, scale);
  }
  fclose(input);
  create_kdtrees();
}
