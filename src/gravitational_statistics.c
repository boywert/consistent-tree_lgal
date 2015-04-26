#undef _POSIX_SOURCE
#include <stdio.h>
#define _POSIX_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <inttypes.h>
#include "gravitational_statistics.h"
#include "gravitational_consistency_subs.h"
#include "gravitational_consistency.h"
#include "universe_time.h"
#include "check_syscalls.h"

#define METRIC_BPDEX 4

extern FILE *logfile;
extern float max_mvir, min_mvir, box_size;
int64_t counts_min_mass_bin=0, counts_max_mass_bin=0, counts_mass_bins=0;
int64_t *counts[NUM_STATS] = {0};

int64_t metric_min_mass_bin=0, metric_max_mass_bin=0, metric_mass_bins=0;
float x_and_v_correlation = 0;
struct metric_struct {
  double avg, sd, max, min, corr, vec_corr, var, tol;
  int64_t count, fixed;
};
struct metric_struct *sigma_x=NULL, *sigma_v=NULL, *sigma_vmax=NULL;
struct metric_struct *sigma_x_subs=NULL, *sigma_v_subs=NULL, *sigma_vmax_subs=NULL;

FILE *metric_output = NULL;

void clear_stats(void)
{
  int64_t i, size;
  metric_min_mass_bin = 
    (min_mvir>0) ? ((int64_t)(METRIC_BPDEX*log10f(min_mvir))) : 0;
  metric_max_mass_bin = 
    (max_mvir>0) ? ((int64_t)(METRIC_BPDEX*log10f(max_mvir))) : 0;
  counts_min_mass_bin = metric_min_mass_bin / METRIC_BPDEX;
  counts_max_mass_bin = metric_max_mass_bin / METRIC_BPDEX;
  counts_mass_bins = counts_max_mass_bin - counts_min_mass_bin + 1;
  metric_mass_bins = metric_max_mass_bin - metric_min_mass_bin + 1;
  assert(counts_min_mass_bin <= counts_max_mass_bin);
  assert(metric_min_mass_bin <= metric_max_mass_bin);

  for (i=0; i<NUM_STATS; i++) {
    counts[i] = check_realloc(counts[i], sizeof(int64_t)*counts_mass_bins, 
			      "consistency statistics");
    memset(counts[i], 0, sizeof(int64_t)*counts_mass_bins);
  }

  size = sizeof(struct metric_struct)*metric_mass_bins;
#define ALLOC_METRIC_STRUCT(x) x = (struct metric_struct *) \
    check_realloc(x, size, "allocating metric struct"); \
  memset(x, 0, size);

  ALLOC_METRIC_STRUCT(sigma_x);
  ALLOC_METRIC_STRUCT(sigma_v);
  ALLOC_METRIC_STRUCT(sigma_vmax);
  ALLOC_METRIC_STRUCT(sigma_x_subs);
  ALLOC_METRIC_STRUCT(sigma_v_subs);
  ALLOC_METRIC_STRUCT(sigma_vmax_subs);
}

inline double calc_dr2(struct tree_halo *h1, struct tree_halo *h2) {
  int64_t i;
  double dx, dr2=0;
  for (i=0; i<3; i++) {
    dx = h1->pos[i] - h2->pos[i];
    if (dx>box_size/2.0) dx -= box_size;
    if (dx<-box_size/2.0) dx += box_size;
    dr2+=dx*dx;
  }
  return dr2;
}

inline double calc_dv2(struct tree_halo *h1, struct tree_halo *h2) {
  int64_t i;
  double dv, dv2=0;
  for (i=0; i<3; i++) {
    dv = h1->vel[i] - h2->vel[i];
    dv2+=dv*dv;
  }
  return dv2;
}

inline float calc_dyn_friction(struct tree_halo *prog, struct tree_halo *parent, struct tree_halo *desc) {
  int64_t i;
  float verr, dv, verr2=0, dv2=0, sum=0, norm;
  for (i=0; i<3; i++) {
    verr = prog->vel[i] - desc->vel[i];
    dv = prog->vel[i] - parent->vel[i];
    sum += verr*dv;
    dv2 += dv*dv;
    verr2 += verr*verr;
  }
  norm = sqrtf(verr2*dv2);
  if (norm) return (sum / norm);
  return 0;
}

inline int64_t metric_bin(float mass) {
  int64_t bin;
  if (!(mass>0)) return -1;
  bin = log10f(mass)*METRIC_BPDEX;
  if ((bin < metric_min_mass_bin) || (bin > metric_max_mass_bin)) return -1;
  return (bin-metric_min_mass_bin);
}

void turn_on_full_metric_output(int64_t output_num, float a1, float a2) {
  float dt = fabs(scale_to_years(a1)-scale_to_years(a2)) / 1.0e6;
  char buffer[1024];
  sprintf(buffer, "%s/full_metric_data_%"PRId64".list", OUTBASE, output_num);
  metric_output = check_fopen(buffer, "w");
  fprintf(metric_output, "#Dt: %f Myr\n", dt);
  fprintf(metric_output, "#From scale %f to %f\n", a1, a2);
  fprintf(metric_output, "#ProgId ProgMass ProgVmax DescId DescMass DescVmax DeltaX DeltaV DeltaVmax Vrms/Vmax DynFricCoeff\n");
}

inline void add_stat(struct metric_struct *m, double val)
{
  //Online variance calculation due to Knuth/Welford
  double delta = val-m->avg;
  m->count++;
  m->avg+=delta / ((double)m->count);
  m->var+=delta*(val - m->avg);
  if (m->max < val) m->max = val;
  if (m->min > val || !m->min) m->min = val;
}

void build_metric_stats(struct tree_halo *prog, struct tree_halo *parent, struct tree_halo *desc)
{
  int64_t bin;
  float dr, dv, dfric=0;
  float dvmax, v_ratio;
  
  if (!(prog->mvir > 0) || !(desc->mvir > 0)
      || !(prog->vmax > 0) || !(desc->vmax > 0)) return; //Unphysical halo
  bin = metric_bin(desc->mvir);
  if (bin < 0) return;

  dvmax = log10f(prog->vmax/desc->vmax);
  dr = sqrtf(calc_dr2(prog, desc));
  dv = sqrtf(calc_dv2(prog, desc));

  add_stat(&(sigma_x[bin]), dr);
  add_stat(&(sigma_v[bin]), dv);
  add_stat(&(sigma_vmax[bin]), dvmax);
  sigma_x[bin].corr += dr*dv;
  
  if (parent) {
    add_stat(&(sigma_x_subs[bin]), dr);
    add_stat(&(sigma_v_subs[bin]), dv);
    add_stat(&(sigma_vmax_subs[bin]), dvmax);
    dfric = calc_dyn_friction(prog, parent, desc);
    sigma_x_subs[bin].corr += dr*dv;
    sigma_v_subs[bin].corr += dfric;
  }


  if (metric_output) {
    v_ratio = desc->vrms / desc->vmax;
    fprintf(metric_output, "%"PRId64" %e %f %"PRId64" %e %f %e %e %e %f %.3f\n",
	    prog->id, prog->mvir, prog->vmax, desc->id, desc->mvir, 
	    desc->vmax, dr, dv, dvmax, v_ratio, dfric);
  }
}

void _finish_metric_stats(struct metric_struct *m) {
  if (!m->count) return;
  m->var /= (double)(m->count);
  if (m->var < 0) m->var = 0;
  m->sd = sqrtf(m->var);
  m->tol = m->avg + m->sd;
}

void _print_metric_stats(int64_t output_num, float a1, float a2, char *title, struct metric_struct *x, struct metric_struct *v, struct metric_struct *vmax)
{
  int64_t i;
  float m, nd;
  char buffer[1024];
  float dt = fabs(scale_to_years(a1)-scale_to_years(a2)) / 1.0e6;
  float vol = pow(box_size/h0, 3);
  FILE *output, *output2;
  snprintf(buffer, 1024, "%s/%s_statistics_%"PRId64".list", OUTBASE, title, output_num);
  output = check_fopen(buffer, "w");
  fprintf(output, "#Dt: %f Myr\n", dt);
  fprintf(output, "#From scale %f to %f\n", a1, a2);
  fprintf(output, "#Mass DX_avg DX_sd DX_tol DV_avg DV_sd DV_tol XV_corr DF_corr DVmax_avg DVmax_sd DVmax_tol Count Number_density\n");
  snprintf(buffer, 1024, "%s/%s_maxmin_%"PRId64".list", OUTBASE, title, output_num);
  output2 = check_fopen(buffer, "w");
  fprintf(output2, "#Dt: %f Myr\n", dt);
  fprintf(output2, "#From scale %f to %f\n", a1, a2);
  fprintf(output2, "#Mass DX_min DX_max DV_min DV_max DVmax_min DVmax_max\n");

  for (i=0; i<metric_mass_bins; i++) {
    if (!x[i].count) continue;
    m = ((float)metric_min_mass_bin + i) / ((float)METRIC_BPDEX);
    nd = log10f(((float)(x[i].count*METRIC_BPDEX))/vol);
    fprintf(output, "%.2f %.3e %.3e %.3e %f %f %f %.2f %.2f %.3f %.3f %.3f %"PRId64" %f\n",
	    m, x[i].avg, x[i].sd, x[i].tol,
	    v[i].avg, v[i].sd, v[i].tol,
	    x[i].corr, v[i].corr, vmax[i].avg,
	    vmax[i].sd, vmax[i].tol, x[i].count, nd);
    fprintf(output2, "%.2f %.3e %.3e %.3e %.3e %.3e %.3e\n",
    	    m, x[i].min, x[i].max, v[i].min,
	    v[i].max, vmax[i].min, vmax[i].max);
  }

  fclose(output);
  fclose(output2);
}

void print_metric_stats(int64_t output_num, float a1, float a2) {
  _print_metric_stats(output_num, a1, a2, "metric", sigma_x, sigma_v, sigma_vmax);
  _print_metric_stats(output_num, a1, a2, "metric_subs", sigma_x_subs, sigma_v_subs, sigma_vmax_subs);
}

void _propagate_statistics(struct metric_struct *m, int64_t tolcomp, int64_t min_stats, int64_t dir) {
  int64_t count = 0;
  int64_t i = (dir == 1) ? 0 : metric_mass_bins-1; //, j;
  double delta;
  struct metric_struct n = {0};
  for (; (i>=0) && (i<metric_mass_bins); i+=dir) {
    count += m[i].count;
    if (!count) continue;
    delta = n.avg - m[i].avg;
    n.var = (n.count*n.var + m[i].count*m[i].var)/(double)count 
      + (n.count*m[i].count)*delta*delta / (double)(count*count);
    n.avg = (n.avg*n.count + m[i].count*m[i].avg) / (double) count;
    n.count = count;
    if (count > min_stats && n.var > 0) break;
  }
  if ((count > min_stats) && (m[i].count > min_stats)) i-=dir;
  if (n.var < 0) n.var = 0;
  n.sd = sqrt(n.var);
  /*if (tolcomp)
    for (j=i; (j>=0) && (j<metric_mass_bins); j-=dir) {
      if (n.sd < m[j].sd) n.sd = m[j].sd;
      if (n.avg < m[j].avg) n.avg = m[j].avg;
      }*/

  
  for (; (i>=0) && (i<metric_mass_bins); i-=dir) {
    m[i].avg = n.avg;
    m[i].sd = n.sd;
    m[i].var = n.var;
    m[i].tol = m[i].sd;
    if (tolcomp) m[i].tol += m[i].avg;
    m[i].fixed = 1;
  }
}

void propagate_statistics(struct metric_struct *m, int64_t tolcomp, int64_t min_stats) {
  _propagate_statistics(m, tolcomp, min_stats, -1);
  _propagate_statistics(m, tolcomp, min_stats, 1);
}

void finish_metric_stats(int64_t output_num, float a1, float a2)
{
  int64_t i;
  if (metric_output) {
    fclose(metric_output);
    metric_output = NULL;
  }
  for (i=0; i<metric_mass_bins; i++) {
    _finish_metric_stats(&(sigma_x[i]));
    _finish_metric_stats(&(sigma_v[i]));
    _finish_metric_stats(&(sigma_vmax[i]));
    _finish_metric_stats(&(sigma_x_subs[i]));
    _finish_metric_stats(&(sigma_v_subs[i]));
    _finish_metric_stats(&(sigma_vmax_subs[i]));
    if (sigma_v_subs[i].count)
      sigma_v_subs[i].corr /= (double) sigma_v_subs[i].count;
    sigma_vmax[i].tol = sigma_vmax[i].sd;
    if (!sigma_vmax[i].tol) sigma_vmax[i].tol = 0.1;
    if (sigma_x[i].count && (sigma_x[i].sd*sigma_v[i].sd)>0) {
      sigma_x[i].corr /= (double) sigma_x[i].count;
      sigma_x[i].corr -= sigma_x[i].avg*sigma_v[i].avg;
      sigma_x[i].corr /= sigma_x[i].sd*sigma_v[i].sd;
    } else {
      sigma_x[i].corr = 0;
    }
    if (sigma_x_subs[i].count && (sigma_x_subs[i].sd*sigma_v_subs[i].sd)>0) {
      sigma_x_subs[i].corr /= (double) sigma_x_subs[i].count;
      sigma_x_subs[i].corr -= sigma_x_subs[i].avg*sigma_v_subs[i].avg;
      sigma_x_subs[i].corr /= sigma_x_subs[i].sd*sigma_v_subs[i].sd;
    } else {
      sigma_x_subs[i].corr = 0;
    }
  }

  //Propagate error calculations for bins with low counts
  propagate_statistics(sigma_x, 1, 20);
  propagate_statistics(sigma_v, 1, 20);
  propagate_statistics(sigma_vmax, 0, 100);
  for (i=1; i<metric_mass_bins; i++) {
    if (sigma_x[i].fixed || ((sigma_x[i].sd) > 0 && (sigma_x[i].count > 1))) continue;
    sigma_x[i] = sigma_x[i-1];
    if (!sigma_vmax[i].fixed) sigma_vmax[i] = sigma_vmax[i-1];
    sigma_v[i] = sigma_v[i-1];
    sigma_x[i].count = sigma_v[i].count = sigma_vmax[i].count = 1;
  }
  for (i=metric_mass_bins-2; i>=0; i--) {
    if (sigma_x[i].fixed || ((sigma_x[i].sd) > 0 && (sigma_x[i].count > 1))) continue;
    sigma_x[i] = sigma_x[i+1];
    if (!sigma_vmax[i].fixed) sigma_vmax[i] = sigma_vmax[i+1];
    sigma_v[i] = sigma_v[i+1];
    sigma_x[i].count = sigma_v[i].count = sigma_vmax[i].count = 1;
  }

  print_metric_stats(output_num, a1, a2);
}

inline void sigma_x_v_vmax(struct tree_halo *h, float *s_x, float *s_v,
			   float *vmax_avg, float *s_vmax, float a) {
  if (h->mvir < 0) return;
  if (metric_mass_bins == 1) {
    *s_x = sigma_x[0].tol;
    *s_v = sigma_v[0].tol;
    *s_vmax = sigma_vmax[0].tol;
    *vmax_avg = sigma_vmax[0].avg;
    return;
  }
  float logm = log10f(h->mvir)-1.0/(float)(METRIC_BPDEX*2);
  float f = logm*METRIC_BPDEX;
  int64_t bin = f;
  f -= bin;
  bin -= metric_min_mass_bin;
  if (bin < 0) { bin = 0; f=0; }
  if (bin > (metric_mass_bins-2)) { bin = metric_mass_bins-2; f=1; }
  *s_x = sigma_x[bin].tol + f*(sigma_x[bin+1].tol-sigma_x[bin].tol);
  *s_v = sigma_v[bin].tol + f*(sigma_v[bin+1].tol-sigma_v[bin].tol);
  *s_vmax = sigma_vmax[bin].tol + f*(sigma_vmax[bin+1].tol-sigma_vmax[bin].tol);
  *vmax_avg = sigma_vmax[bin].avg + f*(sigma_vmax[bin+1].avg-sigma_vmax[bin].avg);
  if (!(*s_x > 0 && *s_v > 0 && *s_vmax > 0)) {
    fprintf(stderr, "Error: too few halos at scale factor %f to calculate consistency metric.\n"
	    "Please remove this and all earlier timesteps from the scale file and rerun.\n(%s)\n", a, SCALEFILE);
    exit(1);
  }
}


void print_stats(int64_t out_num) {
  FILE *output;
  char buffer[1024];
  int64_t i,j;
  snprintf(buffer, 1024, "%s/statistics_%"PRId64".list", OUTBASE, out_num);
  output = check_fopen(buffer, "w");
  fprintf(output, "#T: Total number of halos.\n");
  fprintf(output, "#G: Good halos (with consistent descendants).\n");
  fprintf(output, "#D: Dead halos (without consistent descendants).\n");
  fprintf(output, "#ND: Halos initially without descendants in particle merger tree.\n");
  fprintf(output, "#DNF: Halos whose descendants were deleted.\n");
  fprintf(output, "#PGI: Halos whose descendants are gravitationally inconsistent (binned by progenitor mass).\n");
  fprintf(output, "#DGI: Halos whose descendants are gravitationally inconsistent (binned by descendant mass).\n");
  fprintf(output, "#IT: Halos whose descendants are inconsistent with applicable tidal forces.\n");
  fprintf(output, "#PGR: Halos whose descendants were found by gravitational consistency (binned by progenitor mass).\n");
  fprintf(output, "#DGR: Halos whose descendants were found by gravitational consistency (binned by descendant mass).\n");
  fprintf(output, "#TR: Halos whose descendants were identified by tidal forces.\n");
  fprintf(output, "#P: Phantom halos.\n");
  fprintf(output, "#NP: Phantom halos created at this timestep.\n");
  fprintf(output, "#RND: Halos remaining without descendants (not including phantoms).\n");
  fprintf(output, "#MT: Halos which were likely too massive to tidally merge.\n");
  fprintf(output, "#SMMP: MMP links where the mass ratio was greater than %f.\n", MIN_MMP_MASS_RATIO);
  fprintf(output, "#TMP: Tracks which were truncated due to the phantom fraction being higher than %f.\n", MAX_PHANTOM_FRACTION);
  fprintf(output, "#Mass_bin T G D ND DNF PGI DGI IT PGR DGR TR P NP RND MT SMMP TMP\n");

  for (i=0; i<counts_mass_bins; i++) {
    fprintf(output, "%-2"PRId64" ", counts_min_mass_bin+i);
    for (j=0; j<NUM_STATS; j++)
      fprintf(output, "%7"PRId64" ", counts[j][i]);
    fprintf(output, "\n");
  }
  fclose(output);
}

inline void add_to_mass_bin(int64_t stat, float mass) {
  int64_t bin;
  if (!(mass>0)) return;
  bin = log10f(mass);
  if (bin < counts_min_mass_bin || bin > counts_max_mass_bin) return;
  bin -= counts_min_mass_bin;
#pragma omp atomic
  counts[stat][bin]++;
}

void count_good_halos(struct halo_stash *s) {
  int64_t i;
  for (i=0; i<s->num_halos; i++) {
    if (s->halos[i].descid>=0) add_to_mass_bin(STATS_GOOD, s->halos[i].mvir);
    add_to_mass_bin(STATS_TOTAL, s->halos[i].mvir);
  }
}



void log_too_many_phantoms(float a1, float a2, struct tree_halo *prog) {
  add_to_mass_bin(STATS_TOO_MANY_PHANT, prog->mvir);
  fprintf(logfile, "Removed track with too many phantoms (Phantom count: %"PRId64"; tracking count %"PRId64"): (%f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f into halo %f %"PRId64").\n",
	  prog->num_mmp_phantoms, prog->tracked+1, a1, prog->id, prog->mvir,
	  prog->vmax, prog->np, prog->pos[0], prog->pos[1], prog->pos[2], 
	  prog->vel[0], prog->vel[1], prog->vel[2],
	  a2, prog->descid);
}

void log_phantom_halo(float a1, float a2, struct tree_halo *prog, struct tree_halo *desc) {
  add_to_mass_bin(STATS_PHANTOM, prog->mvir);
  if (prog->phantom == 1) add_to_mass_bin(STATS_NEW_PHANTOM, prog->mvir);
  fprintf(logfile, "Added phantom halo (Phantom count: %"PRId64"): (%f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f into halo %f %"PRId64").\n",
	  prog->phantom, a1, prog->id, prog->mvir, prog->vmax, 
	  prog->np, prog->pos[0], prog->pos[1], prog->pos[2], 
	  prog->vel[0], prog->vel[1], prog->vel[2],
	  a2, desc->id);
}

void log_dead_halo(float a, struct tree_halo *halo) {
  add_to_mass_bin(STATS_DEAD, halo->mvir);
  if (!halo->phantom) add_to_mass_bin(STATS_REM_NODESC, halo->mvir);
  fprintf(logfile, "No descendant found, killing halo. (%f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f).\n",
	  a, halo->id, halo->mvir, halo->vmax,
	  halo->np, halo->pos[0], halo->pos[1],
	  halo->pos[2], halo->vel[0], halo->vel[1],
	  halo->vel[2]);
}

void log_tidal_repair(float a1, float a2, struct tree_halo *halo, struct tree_halo *target)
{
  add_to_mass_bin(STATS_TIDAL_REPAIR, halo->mvir);
  fprintf(logfile, "Tidal disruption conceivable (Force: %f): (%f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f into %f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f) (Desc scale, id: %f %"PRId64").\n",
	      halo->tidal_force, a1, halo->id, halo->mvir, halo->vmax,
	      halo->np, halo->pos[0], 
	      halo->pos[1], halo->pos[2], halo->vel[0], halo->vel[1],
	      halo->vel[2], a1, target->id, target->mvir,
	      target->vmax, target->np,
	      target->pos[0], target->pos[1], target->pos[2], 
	      target->vel[0], target->vel[1], target->vel[2], a2, halo->descid);
}

void log_grav_repair(float a1, float a2, float min_metric, struct tree_halo *prog, struct tree_halo *desc) {
  add_to_mass_bin(STATS_PROG_GRAV_REPAIR, prog->mvir);
  add_to_mass_bin(STATS_DESC_GRAV_REPAIR, desc->mvir);
  fprintf(logfile, "Found match for MMP link (Metric: %f): (%f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f -> %f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f).\n",
	    min_metric, a1, prog->id, prog->mvir, prog->vmax,
	    prog->np, prog->pos[0],
	    prog->pos[1], prog->pos[2], prog->vel[0], prog->vel[1],
	    prog->vel[2], a2, desc->id, desc->mvir,
	    desc->vmax, desc->np,
	    desc->pos[0], desc->pos[1], desc->pos[2],
	    desc->vel[0], desc->vel[1], desc->vel[2]);
}

void log_desc_not_found(float a1, float a2, struct tree_halo *prog) {
  add_to_mass_bin(STATS_DESC_NOT_FOUND, prog->mvir);
  fprintf(logfile, "Severed link: could not find descendant halo. (%f %"PRId64" %e !-> %f %"PRId64")\n",
	  a1, prog->id, prog->mvir, a2, prog->descid);
}

void log_no_desc(float a, struct tree_halo *halo) {
  add_to_mass_bin(STATS_NO_DESC, halo->mvir);
  fprintf(logfile, "Halo has no descendant in particle merger trees: (%f %"PRId64" %e)\n",
	  a, halo->id, halo->mvir);
}

void log_not_enough_tidal(float a1, float a2, struct tree_halo *prog, struct tree_halo *desc) {
  add_to_mass_bin(STATS_INSUFFICIENT_TIDAL, prog->mvir);
  fprintf(logfile, "Severed non-MMP link: Insufficient tidal force (%e).  (%f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f !-> %f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f).\n",
	  prog->tidal_force, a1, prog->id, prog->mvir, prog->vmax, 
	  prog->np, prog->pos[0], prog->pos[1], prog->pos[2],
	  prog->vel[0], prog->vel[1], prog->vel[2],
	  a2, desc->id, desc->mvir, desc->vmax, 
	  desc->np, desc->pos[0], desc->pos[1],
	  desc->pos[2], desc->vel[0], desc->vel[1],
	  desc->vel[2]);
}

void log_break_massive_tidal(float a1, float a2, struct tree_halo *prog, struct tree_halo *desc) {
  add_to_mass_bin(STATS_MASSIVE_TIDAL, prog->mvir);
  fprintf(logfile, "Severed non-MMP link: Likely that halo should still be tracked (num particles: %"PRId64").  (%f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f !-> %f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f).\n",
	  prog->np, a1, prog->id, prog->mvir, prog->vmax, 
	  prog->np, prog->pos[0], prog->pos[1], prog->pos[2],
	  prog->vel[0], prog->vel[1], prog->vel[2],
	  a2, desc->id, desc->mvir, desc->vmax, 
	  desc->np, desc->pos[0], desc->pos[1],
	  desc->pos[2], desc->vel[0], desc->vel[1],
	  desc->vel[2]);
}

void log_spurious_mmp_link_broken(float a1, float a2, struct tree_halo *prog, struct tree_halo *desc) {
  add_to_mass_bin(STATS_SPURIOUS_MMP, prog->mvir);
  fprintf(logfile, "Severed MMP link: Mass ratio greater than %f.  (%f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f !-> %f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f).\n",
	  MIN_MMP_MASS_RATIO, a1, prog->id, prog->mvir, prog->vmax, 
	  prog->np, prog->pos[0], prog->pos[1], prog->pos[2],
	  prog->vel[0], prog->vel[1], prog->vel[2],
	  a2, desc->id, desc->mvir, desc->vmax, 
	  desc->np, desc->pos[0], desc->pos[1],
	  desc->pos[2], desc->vel[0], desc->vel[1],
	  desc->vel[2]);
}

void log_grav_inconsistency(float a1, float a2, float metric_val, struct tree_halo *prog, struct tree_halo *desc) {
  add_to_mass_bin(STATS_PROG_GRAV_INCONSISTENCY, prog->mvir);
  add_to_mass_bin(STATS_DESC_GRAV_INCONSISTENCY, desc->mvir);
  fprintf(logfile, "Severed MMP link: Possibly inconsistent with law of gravity (Metric: %f).  (%f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f !-> %f %"PRId64" %e %f %"PRId64" %f %f %f %f %f %f).\n",
	  metric_val, a1, prog->id, prog->mvir, prog->vmax, prog->np,
	  prog->pos[0], prog->pos[1], prog->pos[2], prog->vel[0],
	  prog->vel[1], prog->vel[2], a2, desc->id,
	  desc->mvir, desc->vmax, desc->np,
	  desc->pos[0], desc->pos[1], desc->pos[2], 
	  desc->vel[0], desc->vel[1], desc->vel[2]);
}
