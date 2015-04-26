#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include <strings.h>
#include <sys/stat.h>
#include <time.h>
#include "gravitational_consistency.h"
#include "masses.h"
#include "check_syscalls.h"
#include "stringparse.h"

extern float min_mvir, max_mvir;
extern float box_size;
FILE *timing_file = NULL;
int64_t timing_start = 0;

void print_timing(char *prog, char *info) {
  char buffer[1024];
  if (!timing_file) {
    snprintf(buffer, 1024, "%s/timing.log", OUTBASE);
    timing_file = check_fopen(buffer, "a");
    timing_start = time(NULL);
    time_t t = timing_start;
    fprintf(timing_file, "%s started at %s", prog, ctime(&t));
    fflush(timing_file);
  }
  if (!info) return;
  fprintf(timing_file, "[%6"PRId64"s] %s\n", (int64_t)time(NULL)-timing_start, info);
  fflush(timing_file);
}

void close_timing_log(void) {
  fprintf(timing_file, "\n");
  fclose(timing_file);
}


void gzip_file(char *filename) {
  char buffer[1024];
  snprintf(buffer, 1024, "gzip -f \"%s\" &", filename);
  if (!LIMITED_MEMORY)
    system(buffer);
}

void print_halo(FILE *o, struct tree_halo *th) {
  int64_t mmp, i;
  int64_t flags;
  if (!th) {
    if (!strcasecmp(INPUT_FORMAT, "BINARY")) return;
    fprintf(o, "#ID DescID Mvir Vmax Vrms Rvir Rs Np X Y Z VX VY VZ Jx Jy Jz spin %s Phantom MMP Suspicious? PID UPID Tracked Tracked_Single_MMP Num_MMP_Phantoms Original_ID Last_mm\n", EXTRA_PARAM_LABELS);
     fprintf(o, "#Units: Masses in Msun / h\n"
            "#Units: Positions in Mpc / h (comoving)\n"
            "#Units: Velocities in km / s (physical)\n"
            "#Units: Angular Momenta in (Msun/h) * (Mpc/h) * km/s (physical)\n"
            "#Units: Radii in kpc / h (comoving)\n");
    return;
  }
  if (!strcasecmp(INPUT_FORMAT, "BINARY")) {
    fwrite(th, sizeof(struct tree_halo), 1, o);
    return;
 }
  mmp = (th->flags & MMP_FLAG) ? 1 : 0;
  flags = 0;
  if (th->flags & SUSPICIOUS_LINK_FLAG) flags = 1;
  if (th->flags & UNPHYSICAL_LINK_FLAG) flags = 2;
  if (th->flags & MERGER_FLUCTUATION_FLAG) flags = 3;
  flags += (th->flags & (TOO_MANY_PHANT_FLAG | SHORT_TRACK_FLAG | MISTRACKED_SUB_FLAG));
  fprintf(o, "%"PRId64" %"PRId64" %.3e %.2f %.2f %.3f %.3f %"PRId64" %.5f %.5f %.5f %.2f %.2f %.2f %.3e %.3e %.3e %.5f",
	  th->id, th->descid, th->mvir, th->vmax, th->vrms, th->rvir, th->rs,
	  th->np, th->pos[0], th->pos[1], th->pos[2], th->vel[0], th->vel[1],
	  th->vel[2], th->J[0], th->J[1], th->J[2], th->spin);
  for (i=0; i<EXTRA_PARAMS; i++) {
    if (th->extra_params[i] == ((double)((int64_t)th->extra_params[i])))
      fprintf(o, " %"PRId64, (int64_t)th->extra_params[i]);
    else
      fprintf(o, " %e", th->extra_params[i]);
  }
  fprintf(o, " %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %f\n",
	  th->phantom, mmp, flags, th->pid, th->upid, th->tracked,
	  th->tracked_single_mmp, th->num_mmp_phantoms, th->orig_id, th->last_mm);
}

void common_print_halos(char *filename, struct tree_halo *h, int64_t num_halos, int64_t gzip_halos)
{
  int64_t i;
  FILE *o;

  o = check_fopen(filename, "w");
  if (!strcasecmp(INPUT_FORMAT, "BINARY")) {
    fwrite(h, sizeof(struct tree_halo), num_halos, o);
  } else {
    print_halo(o, NULL);
    for (i=0; i<num_halos; i++) print_halo(o, &(h[i]));
  }
  fclose(o);
  if (gzip_halos) gzip_file(filename);
}

int really_beyond_mmp_ratio(struct tree_halo h1, struct tree_halo h2) {
  float m1, m2, v1, v2;
  if (h1.mvir > h2.mvir) { m1 = h1.mvir; m2 = h2.mvir; }
  else { m1 = h2.mvir; m2 = h1.mvir; }
  if (m1*0.3*MIN_MMP_MASS_RATIO > m2) return 1;

  if (h1.vmax > h2.vmax) { v1 = h1.vmax; v2 = h2.vmax; }
  else { v2 = h1.vmax; v1 = h2.vmax; }
  if (v1*0.5*MIN_MMP_VMAX_RATIO > v2) return 1;

  return 0;
}

int beyond_mmp_ratio(struct tree_halo h1, struct tree_halo h2) {
  float m1, m2, v1, v2;
  if (h1.mvir > h2.mvir) { m1 = h1.mvir; m2 = h2.mvir; }
  else { m1 = h2.mvir; m2 = h1.mvir; }
  if (m1*MIN_MMP_MASS_RATIO > m2) return 1;

  if (h1.vmax > h2.vmax) { v1 = h1.vmax; v2 = h2.vmax; }
  else { v2 = h1.vmax; v1 = h2.vmax; }
  if (v1*MIN_MMP_VMAX_RATIO > v2) return 1;

  return 0;
}

void build_id_conv_list(struct halo_stash *h)
{
  int64_t n, ids;
  if (h->num_halos < 1) return;

  h->min_id = h->max_id = h->halos[0].id;
  for (n=0; n<h->num_halos; n++) {
    if (h->max_id < h->halos[n].id) h->max_id = h->halos[n].id;
    if (h->min_id > h->halos[n].id) h->min_id = h->halos[n].id;
  }

  ids = (1 + h->max_id - h->min_id);
  h->id_conv = (int64_t *)check_realloc(h->id_conv, sizeof(int64_t)*ids, "halo index conversion");
  for (n=0; n<ids; n++) h->id_conv[n]=-1;

  for (n=0; n<h->num_halos; n++)
    h->id_conv[h->halos[n].id - h->min_id] = n;

  assert(h->min_id >= 0); //Otherwise, logic problems occur in rest of code.
}


void add_halo(struct halo_stash *h, struct tree_halo d)
{
  if (((h->num_halos)%1000)==0) {
    int64_t array_size = h->num_halos+1000;
    h->halos = (struct tree_halo *)
      check_realloc(h->halos,sizeof(struct tree_halo)*array_size, "Adding new halo.");
  }

  h->halos[h->num_halos] = d;
  h->num_halos++;
}

void load_halos(char *filename, struct halo_stash *h, float scale, int dead)
{
  FILE *input;
  char buffer[1024];
  struct tree_halo d = {0};
  int64_t n;
  int64_t mmp, flags, phantom;
  float delta_mvir = delta_vir(scale);
  float mean_density = 2.77519737e11*Om; //(Msun/h) / (comoving Mpc/h)^3
  float vir_density = delta_mvir*mean_density;
  int64_t expected_inputs = 28+EXTRA_PARAMS;
  int64_t regular_inputs = 18;
  SHORT_PARSETYPE;
#define NUM_INPUTS (28+MAX_EXTRA_PARAMS)
  enum short_parsetype stypes[NUM_INPUTS] = 
    { D64, D64, F, F, F, F, F, D64, F, F, F, F, F, F, F, F, F, F, D64, D64, D64, D64, D64, D64, D64, D64, D64, F };
  enum parsetype types[NUM_INPUTS];
  void *data[NUM_INPUTS] = {&(d.id),
			    &(d.descid), &(d.mvir), &(d.vmax), &(d.vrms),
			    &(d.rvir), &(d.rs), &(d.np), &(d.pos[0]),
			    &(d.pos[1]), &(d.pos[2]), &(d.vel[0]), 
			    &(d.vel[1]), &(d.vel[2]), &(d.J[0]), &(d.J[1]),
			    &(d.J[2]), &(d.spin),
			    &(phantom), &(mmp), &(flags), &(d.pid), &(d.upid),
			    &(d.tracked), &(d.tracked_single_mmp), 
			    &(d.num_mmp_phantoms), &(d.orig_id), &(d.last_mm)};

  int64_t read_ascii = 1;
  if (!strcasecmp(INPUT_FORMAT, "BINARY")) read_ascii = 0;

  gen_ff_cache();
  d.num_prog = d.flags = d.phantom = 0;
  d.mmp_id = -1;
  for (n=0; n<NUM_INPUTS; n++) types[n] = stypes[n];
  for (n=expected_inputs-1; n>=regular_inputs+EXTRA_PARAMS; n--) {
    types[n] = types[n-EXTRA_PARAMS];
    data[n] = data[n-EXTRA_PARAMS];
  }
  for (n=regular_inputs; n<regular_inputs+EXTRA_PARAMS; n++) {
    types[n] = PARSE_FLOAT64;
    data[n] = &(d.extra_params[n-regular_inputs]);
  }
  
  int64_t old_halos = h->num_halos;
  input = check_fopen(filename, "r");
  if (read_ascii) {
    while (fgets(buffer, 1024, input)) {
      if (buffer[0] == '#') continue;
      mmp = flags = 0;
      n = stringparse(buffer, data, (enum parsetype *)types, expected_inputs);
    
      if (n < regular_inputs-4) continue;
      if (n < regular_inputs) { d.J[0] = d.J[1] = d.J[2] = d.spin = 0; }
      if (n < regular_inputs+EXTRA_PARAMS+1) d.phantom = 0;
      else d.phantom = phantom;
      if (n < regular_inputs+EXTRA_PARAMS+4) d.pid = -1;
      if (n < regular_inputs+EXTRA_PARAMS+5) d.upid = -1;
      if (n < regular_inputs+EXTRA_PARAMS+6) d.tracked = 0;
      if (n < regular_inputs+EXTRA_PARAMS+7) d.tracked_single_mmp = 0;
      if (n < regular_inputs+EXTRA_PARAMS+8) d.num_mmp_phantoms = 0;
      if (n < regular_inputs+EXTRA_PARAMS+9) d.orig_id = -1;
      if (n < regular_inputs+EXTRA_PARAMS+10) d.last_mm = 0;

      d.mvir = fabs(d.mvir);
      if (!(d.mvir>0)) continue;
      if (!(d.rvir > 0))
	d.rvir = cbrt(d.mvir / (4.0*M_PI*vir_density/3.0)) * 1000.0;
      if (!(d.rs > 0)) d.rs = d.rvir / concentration(d.mvir, scale);

      d.flags = 0;
      if (mmp) d.flags |= MMP_FLAG;
      if ((flags & 3)==1) d.flags |= SUSPICIOUS_LINK_FLAG;
      if ((flags & 3)==2) d.flags |= UNPHYSICAL_LINK_FLAG;
      if (dead) d.flags |= DEAD_HALO_FLAG;
      add_halo(h, d);
    }
  }
  else {
    struct stat s;
    assert(!fstat(fileno(input), &s));
    int64_t num_inputs = s.st_size / (sizeof(struct tree_halo));
    if (num_inputs * sizeof(struct tree_halo) != s.st_size) {
      fprintf(stderr, "Error: file format in %s not recognized as %s!\n", filename, INPUT_FORMAT);
      exit(1);
    }
    int64_t alloc_size = num_inputs+1000;    
    h->halos = (struct tree_halo *)
      check_realloc(h->halos,sizeof(struct tree_halo)*(old_halos+alloc_size), "Adding new halos.");
    fread(h->halos+old_halos, sizeof(struct tree_halo), num_inputs, input);
    h->num_halos += num_inputs;
  }
  fclose(input);

  for (int64_t i=old_halos; i<h->num_halos; i++) {
    mmp = h->halos[i].flags & MMP_FLAG;
    flags = 0;
    if (h->halos[i].flags & SUSPICIOUS_LINK_FLAG) flags |= 1;
    if (h->halos[i].flags & UNPHYSICAL_LINK_FLAG) flags |= 2;
    h->halos[i].flags = 0;
    if (mmp) h->halos[i].flags |= MMP_FLAG;
    if ((flags & 3)==1) h->halos[i].flags |= SUSPICIOUS_LINK_FLAG;
    if ((flags & 3)==2) h->halos[i].flags |= UNPHYSICAL_LINK_FLAG;
    if (dead) h->halos[i].flags |= DEAD_HALO_FLAG;
    if (min_mvir==0) min_mvir = h->halos[i].mvir;
    if (h->halos[i].mvir < min_mvir) min_mvir = h->halos[i].mvir;
    if (h->halos[i].mvir > max_mvir) max_mvir = h->halos[i].mvir;
    for (n=0; n<3; n++) if (h->halos[i].pos[n] > box_size) box_size = h->halos[i].pos[n];
  }

  if (!BOX_WIDTH)
    box_size = (int) (box_size + 0.5);
  else
    box_size = BOX_WIDTH;
}


void read_outputs(float **output_scales, int64_t **outputs, int64_t *num_outputs) {
  FILE *input;
  int64_t n;
  char buffer[1024];
  
  input = check_fopen(SCALEFILE, "r");
  while (fgets(buffer, 1024, input)) {
    if (((*num_outputs)%1000 == 0)) {
      *output_scales = (float *)check_realloc(*output_scales, sizeof(float)*((*num_outputs)+1000), "output scales");
      *outputs = (int64_t *)check_realloc(*outputs, sizeof(int64_t)*((*num_outputs)+1000), "output numbers");
    }
    n = sscanf(buffer, "%"SCNd64" %f", &((*outputs)[*num_outputs]), &((*output_scales)[*num_outputs]));
    if (n<2) continue;
    if (num_outputs[0]) {
      float scale_now = output_scales[0][num_outputs[0]];
      float last_scale = output_scales[0][num_outputs[0]-1];
      if (scale_now == last_scale) continue;
      if (scale_now < last_scale) {
	fprintf(stderr, "[Error] Scale file (%s) not sorted in ascending order!\n", SCALEFILE);
	exit(1);
      }
    }
    *num_outputs = (*num_outputs) + 1;
  }
  fclose(input);

  *output_scales = (float *)check_realloc(*output_scales, sizeof(float)*((*num_outputs)), "output scales");
  *outputs = (int64_t *)check_realloc(*outputs, sizeof(int64_t)*((*num_outputs)), "output numbers");
}
