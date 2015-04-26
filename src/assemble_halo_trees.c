#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <sys/resource.h>
#include "halo_io.h"
#include "cached_io.h"
#include "litehash.h"
#include "grav_config.h"
#include "gravitational_consistency.h"
#include "gravitational_consistency_vars.h"
#include "stringparse.h"
#include "assemble_halo_trees.h"
#include "check_syscalls.h"
#include "version.h"

FILE **tree_outputs = NULL;
float *output_scales=NULL;
int64_t *outputs=NULL, num_outputs=0;
int64_t *num_trees = NULL;
struct merger_halo *extra_halos = NULL;
FILE **tree_inputs = NULL;
struct merger_halo *halos = NULL;
struct litehash *lh = NULL;
int64_t num_halos=0, num_halos_output=0, num_trees_output = 0;
float box_size, min_mvir, max_mvir;
FILE *locations = NULL;

int main(int argc, char **argv)
{
  int64_t i, j, k;
  int64_t tree_header_offset;
  char buffer[1024];
  //struct cached_io *input;
  FILE *input;
  FILE *output;
  struct rlimit rlp;

  getrlimit(RLIMIT_NOFILE, &rlp);
  rlp.rlim_cur = rlp.rlim_max;
  setrlimit(RLIMIT_NOFILE, &rlp);

  if (argc > 1) grav_config(argv[1], 0);
  else { 
    fprintf(stderr, "Consistent Trees, Version %s\n", TREE_VERSION);
    fprintf(stderr, "%s.  See the LICENSE file for redistribution details.\n", TREE_COPYRIGHT);
    fprintf(stderr, "Usage: %s options.cfg\n", argv[0]); exit(1); }
  read_outputs(&output_scales, &outputs, &num_outputs);
  print_timing(_AHT, NULL);

  tree_header_offset = create_headers();
  /*tree_inputs = check_realloc(NULL, sizeof(struct cached_io *)*num_outputs,
    "Allocating tree inputs.");*/
  tree_inputs = check_realloc(NULL, sizeof(FILE *)*num_outputs,
			      "Allocating tree inputs.");
  extra_halos = check_realloc(NULL, sizeof(struct merger_halo)*num_outputs,
			      "Allocating extra halo slots.");
  for (i=0; i<num_outputs; i++) extra_halos[i].id = -1;
  for (i=0; i<num_outputs; i++) {
    snprintf(buffer, 1024, "%s/really_consistent_%"PRId64".list", 
	     OUTBASE, outputs[i]);
    //tree_inputs[i] = cfopen(buffer, 2000000);
    tree_inputs[i] = check_fopen(buffer, "r");
  }
  snprintf(buffer, 1024, "%s/locations.dat", TREE_OUTBASE);
  locations = check_fopen(buffer, "w");
  fprintf(locations, "#TreeRootID FileID Offset Filename\n");

  input = tree_inputs[num_outputs-1];

  struct merger_halo halo;
  int64_t last_timing_output = 0;
  while (read_halo_from_file(&halo, input, num_outputs-1)) {
    halo.scale = output_scales[num_outputs-1];
    halo.desc_scale = 0;

    int64_t id = halo_pos_to_tree(&halo);
    num_trees[id]++;

    num_halos = 0;
    if (lh) { free_litehash2(lh); }
    lh = new_litehash(8);
    halos = add_to_array(halos, &num_halos, sizeof(struct merger_halo), &halo);

    build_tree(halo.id, num_outputs-2);
    print_tree_halos(id);
    num_halos_output += num_halos;
    num_trees_output++;
    if (num_halos_output >= last_timing_output+1000000) {
      timed_output(_AHT, "Wrote %"PRId64" trees (%"PRId64" halos).", 
		   num_trees_output, num_halos_output);
      last_timing_output = 1000000*((int64_t)(num_halos_output/1000000));
    }
  }
  timed_output(_AHT, "Wrote %"PRId64" trees (%"PRId64" halos).",
	       num_trees_output, num_halos_output);

  //Print num_tree counts
  for (i=0; i<BOX_DIVISIONS; i++) {
    for (j=0; j<BOX_DIVISIONS; j++) {
      for (k=0; k<BOX_DIVISIONS; k++) {
	int64_t id = i*BOX_DIVISIONS*BOX_DIVISIONS + j*BOX_DIVISIONS + k;
	output = tree_outputs[id];
	fseek(output, tree_header_offset, SEEK_SET);
	fprintf(output, "%-18"PRId64"\n", num_trees[id]);
	fclose(output);
      }
    }
  }
  fclose(locations);
  timed_output(_AHT, "Successfully finished.");
  close_timing_log();
  return 0;
}

void calc_depthfirst_order(void) {
  int64_t *cur_id = NULL;
  int64_t num_scales = 0;
  double last_scale = -1;
  int64_t depthfirst_order = 0, i, cur_halo = 0, cur_scale = 0;
  for (i=0; i<num_halos; i++) {
    if (halos[i].scale != last_scale) {
      if (!(num_scales % 100))
	cur_id = check_realloc(cur_id, sizeof(int64_t)*(num_scales+100),
			       "Allocating scales");
      cur_id[num_scales] = i;
      num_scales++;
      last_scale = halos[i].scale;
    }
    halos[i].depthfirst_order = -1;
  }

  for (; depthfirst_order < num_halos; depthfirst_order++) {
    if (halos[cur_halo].depthfirst_order == -1)
      halos[cur_halo].depthfirst_order = depthfirst_order;
    else depthfirst_order--;
    //Check progenitor:
    if ((cur_scale < num_scales-1) && (cur_id[cur_scale+1] < num_halos)) {
      double desc_scale = halos[cur_id[cur_scale+1]].desc_scale;
      int64_t desc_id = halos[cur_id[cur_scale+1]].descid;
      struct merger_halo *desc = lh_getval2(lh, &desc_scale, &desc_id);
      if (desc->id == halos[cur_halo].id) {
	cur_scale++;
	cur_halo = cur_id[cur_scale];
	continue;
      }
    }

    //Otherwise, exhausted all progenitors and need to check descendants:
    if ((cur_id[cur_scale]+1) < num_halos) {
      if (halos[cur_id[cur_scale]].scale != halos[cur_id[cur_scale]+1].scale)
	cur_id[cur_scale] = num_halos;
      else cur_id[cur_scale]++;
    }
    else cur_id[cur_scale] = num_halos;

    if (cur_scale > 0) cur_scale--;
    else assert(depthfirst_order == num_halos-1);
    cur_halo = cur_id[cur_scale];
  }

  if (num_halos) free(cur_id);
}

void calc_progenitor_links(void) {
  int64_t i;
  for (i=0; i<num_halos; i++) {
    if ((i+1 < num_halos) && (halos[i+1].descid == halos[i].descid))
      halos[i].next_coprogenitor_df = halos[i+1].depthfirst_order;
    else
      halos[i].next_coprogenitor_df = -1;
    halos[i].last_progenitor_df = halos[i].depthfirst_order;
  }
  for (i=num_halos-1; i>0; i--) {
    struct merger_halo *desc=lh_getval2(lh,&(halos[i].desc_scale),&(halos[i].descid));
    if (desc && desc->last_progenitor_df < halos[i].last_progenitor_df)
      desc->last_progenitor_df = halos[i].last_progenitor_df;
  }
}

void print_tree_halos(int64_t file_id) {
  int64_t i,j,k,offset;
  for (i=0; i<num_halos; i++) {
    lh_setval2(lh, &(halos[i].scale), &(halos[i].id), &(halos[i]));
    halos[i].treeroot_id = halos[0].id;
    halos[i].breadthfirst_order = i;
  }

  calc_depthfirst_order();
  calc_progenitor_links();
  calc_desc_pid_and_num_progs();
  calc_mmp();
  smooth_mass();
  conserve_mass();
  calc_mmp();

  fprintf(tree_outputs[file_id], "#tree %"PRId64"\n", halos[0].id);
  i = file_id/(BOX_DIVISIONS*BOX_DIVISIONS);
  j = (file_id%((int64_t)(BOX_DIVISIONS*BOX_DIVISIONS)))/BOX_DIVISIONS;
  k = file_id%((int64_t)BOX_DIVISIONS);
  offset = ftello(tree_outputs[file_id]);
  fprintf(locations, "%"PRId64" %"PRId64" %"PRId64" tree_%"PRId64"_%"PRId64"_%"PRId64".dat\n",
	  halos[0].id, file_id, offset, i, j, k);
  for (i=0; i<num_halos; i++)
    print_tree_halo(halos + i, tree_outputs[file_id]);
}


void calc_desc_pid_and_num_progs(void) {
  int64_t i;
  for (i=0; i<num_halos; i++) {
    double desc_scale = halos[i].desc_scale;
    int64_t desc_id = halos[i].descid;
    struct merger_halo *desc = lh_getval2(lh, &desc_scale, &desc_id);
    if (!desc) continue;
    halos[i].desc = desc;
    desc->num_prog++;
    halos[i].desc_pid = desc->pid;
  }
}

void calc_mmp(void) {
  int64_t i;
  for (i=0; i<num_halos; i++) {
    struct merger_halo *th = halos + i;
    if (!th->desc) continue;
    if (!th->desc->mmp_halo) th->desc->mmp_halo = th;
    else if (th->desc->mmp_halo->mvir < th->mvir)
      th->desc->mmp_halo = th;
  }

  for (i=0; i<num_halos; i++) {
    if (!halos[i].desc) { halos[i].mmp = 1; continue; }
    halos[i].mmp = (halos[i].desc->mmp_halo->id == halos[i].id) ? 1 : 0;
  }
}

void smooth_mass(void) {
  int64_t i;
  for (i=0; i<num_halos; i++) {
    struct merger_halo *th = halos+i;
    if (!th->desc) continue;
    if (th->desc && th->desc->mmp_halo && (th->desc->mmp_halo->id == th->id)) {
      th->next_mass = th->desc->mvir;
      th->desc->prev_mass = th->mvir;
    }
  }
  for (i=0; i<num_halos; i++) {
    struct merger_halo *th = halos+i;
    if (th->next_mass && th->prev_mass)
      th->mvir = (th->next_mass + th->prev_mass + th->mvir) / 3.0;
  }
}

int double_sort(const void *a, const void *b) {
  double c = *((double *)a);
  double d = *((double *)b);
  if (c<d) return -1;
  if (c>d) return 1;
  return 0;
}

void conserve_mass(void) {
  int64_t i, j, *ids;
  double scale, *scales = lh_keylist(lh);
  qsort(scales, lh->elems, sizeof(double), double_sort);
  for (i=0; i<num_halos; i++) halos[i].incoming_mass = 0;

  for (i=0; i<lh->elems; i++) {
    scale = scales[i];
    struct litehash *idlh = lh_getval(lh, &scale);
    assert(idlh);
    ids = lh_keylist(idlh);
    for (j=0; j<idlh->elems; j++) {
      struct merger_halo *th = lh_getval(idlh, ids+j);
      assert(th);
      if ((th->upid < 0) && th->incoming_mass
	  && (th->incoming_mass > th->mvir)) {
	th->mvir = th->incoming_mass;
      }
      if (!th->desc) continue;
	    
      if (th->pid < 0) th->desc->incoming_mass += th->mvir;
      if ((th->pid < 0) && (th->desc->pid > -1) && (th->desc->id != th->desc->pid)) {
	struct merger_halo *parent = lh_getval2(lh, &(th->desc->scale), &(th->desc->pid));
	if (parent) parent->incoming_mass += th->mvir;
      }
    }
    free(ids);
  }
  free(scales);
}

int64_t create_headers(void) {
  int64_t i, j, k, id, offset = 0;
  char buffer[1024];
  FILE *output;
  int64_t num_output_trees = BOX_DIVISIONS*BOX_DIVISIONS*BOX_DIVISIONS;
  char *epd_prefix = "", *epd_postfix = "", *epd = "";
  tree_outputs = check_realloc(tree_outputs, sizeof(FILE *)*num_output_trees,
			       "Allocating outputs.");

  num_trees = check_realloc(NULL, sizeof(int64_t)*num_output_trees,
			    "Allocating tree counts.");
  memset(num_trees, 0, sizeof(int64_t)*num_output_trees);

  if (EXTRA_PARAMS) {
    epd = EXTRA_PARAM_DESCRIPTIONS;
    if (epd[0] != '#') epd_prefix = "#";
    if (epd[strlen(epd)-1] != '\n') epd_postfix = "\n";
  }

  for (i=0; i<BOX_DIVISIONS; i++) {
    for (j=0; j<BOX_DIVISIONS; j++) {
      for (k=0; k<BOX_DIVISIONS; k++) {
	snprintf(buffer, 1024, "%s/tree_%"PRId64"_%"PRId64"_%"PRId64".dat",
		 TREE_OUTBASE, i, j, k);
	unlink(buffer);
	id = i*BOX_DIVISIONS*BOX_DIVISIONS + j*BOX_DIVISIONS + k;
	output = tree_outputs[id] = check_fopen(buffer, "w");
	fprintf(output,
		"#scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_mvir(9) mvir(10) rvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29) Orig_halo_ID(30) Snap_num(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) %s\n"
		"#Omega_M = %f; Omega_L = %f; h0 = %f\n"
		"#Full box size = %f Mpc/h\n"
		"#Scale: Scale factor of halo.\n"
		"#ID: ID of halo (unique across entire simulation).\n"
		"#Desc_Scale: Scale of descendant halo, if applicable.\n"
		"#Descid: ID of descendant halo, if applicable.\n"
		"#Num_prog: Number of progenitors.\n"
		"#Pid: ID of least massive host halo (-1 if distinct halo).\n"
		"#Upid: ID of most massive host halo (different from Pid when the halo is within two or more larger halos).\n"
		"#Desc_pid: Pid of descendant halo (if applicable).\n"
		"#Phantom: Nonzero for halos interpolated across timesteps.\n"
		"#SAM_Mvir: Halo mass, smoothed across accretion history; always greater than sum of halo masses of contributing progenitors (Msun/h).  Only for use with select semi-analytic models.\n"
		"#Mvir: Halo mass (Msun/h).\n"
		"#Rvir: Halo radius (kpc/h comoving).\n"
		"#Rs: Scale radius (kpc/h comoving).\n"
		"#Vrms: Velocity dispersion (km/s physical).\n"
		"#mmp?: whether the halo is the most massive progenitor or not.\n"
		"#scale_of_last_MM: scale factor of the last major merger (Mass ratio > %g).\n"
		"#Vmax: Maxmimum circular velocity (km/s physical).\n"
		"#X/Y/Z: Halo position (Mpc/h comoving).\n"
		"#VX/VY/VZ: Halo velocity (km/s physical).\n"
		"#JX/JY/JZ: Halo angular momenta ((Msun/h) * (Mpc/h) * km/s (physical)).\n"
		"#Spin: Halo spin parameter.\n"
		"#Breadth_first_ID: breadth-first ordering of halos within a tree.\n"		
		"#Depth_first_ID: depth-first ordering of halos within a tree.\n"
		"#Tree_root_ID: ID of the halo at the last timestep in the tree.\n"
		"#Orig_halo_ID: Original halo ID from halo finder.\n"
		"#Snap_num: Snapshot number from which halo originated.\n"
		"#Next_coprogenitor_depthfirst_ID: Depthfirst ID of next coprogenitor.\n"
		"#Last_progenitor_depthfirst_ID: Depthfirst ID of last progenitor.\n"
		"%s%s%s"
		"#Consistent Trees Version %s\n", EXTRA_PARAM_LABELS, 
		Om, Ol, h0, BOX_WIDTH,
		MAJOR_MERGER, epd_prefix, epd, epd_postfix, TREE_VERSION);
	if (FIX_ROCKSTAR_SPINS > -1) {
	  fprintf(output, "#Includes fix for Rockstar spins & T/|U| (assuming T/|U| = column %"PRId64")\n", FIX_ROCKSTAR_SPINS+34);
	}
	
	offset = ftello(output);
	fprintf(output, "XXXXXXXXXXXXXXXXXX\n"); //For the number of trees
      }
    }
  }
  return offset;
}

int64_t read_halo_from_line(struct merger_halo *halo, char *buffer, int64_t snapnum) {
  SHORT_PARSETYPE;
#define NUM_INPUTS 28+MAX_EXTRA_PARAMS
  enum short_parsetype stypes[NUM_INPUTS] = 
    { D64, D64, F64, F64, F64, F64, F64, D64, F, F, F, F, F, F, F, F, F, F, D64, D64, K, D64, D64, K, K, K, D64, F };
  enum parsetype types[NUM_INPUTS];
  int64_t regular_inputs = 18;
  int64_t expected_inputs = 28+EXTRA_PARAMS;
  int64_t n;
  void *data[NUM_INPUTS] = {&(halo->id),
                            &(halo->descid), &(halo->mvir), &(halo->vmax), 
			    &(halo->vrms), &(halo->rvir), &(halo->rs), 
			    &(halo->np), &(halo->pos[0]), &(halo->pos[1]), 
			    &(halo->pos[2]), &(halo->vel[0]), &(halo->vel[1]), 
			    &(halo->vel[2]), &(halo->J[0]), &(halo->J[1]),
			    &(halo->J[2]), &(halo->spin), &(halo->phantom), 
			    &(halo->mmp), 
			    NULL, &(halo->pid), &(halo->upid), NULL, NULL,
			    NULL, &(halo->orig_id), &(halo->last_mm)};

  memset(halo, 0, sizeof(struct merger_halo));
  for (n=0; n<NUM_INPUTS; n++) types[n] = stypes[n];
  for (n=expected_inputs-1; n>=regular_inputs+EXTRA_PARAMS; n--) {
    types[n] = types[n-EXTRA_PARAMS];
    data[n] = data[n-EXTRA_PARAMS];
  }
  for (n=regular_inputs; n<regular_inputs+EXTRA_PARAMS; n++) {
    types[n] = PARSE_FLOAT64;
    data[n] = &(halo->extra_params[n-regular_inputs]);
  }

  n = stringparse(buffer, data, (enum parsetype *)types, expected_inputs);
  halo->orig_mvir = halo->mvir;
  halo->desc_pid = -1;
  halo->snapnum = snapnum;
  return ((n==expected_inputs) ? 1 : 0);
#undef NUM_INPUTS
}

int64_t read_binary_halo_from_file(struct merger_halo *halo, FILE *input, int64_t snapnum) {
  struct tree_halo hf = {0};
  int64_t n = fread(&hf, sizeof(struct tree_halo), 1, input);
  if (n < 1) return 0;
  memset(halo, 0, sizeof(struct merger_halo));
  halo->id = hf.id;
  halo->descid = hf.descid;
  halo->orig_mvir = halo->mvir = hf.mvir;
  halo->vmax = hf.vmax;
  halo->vrms = hf.vrms;
  halo->rvir = hf.rvir;
  halo->rs = hf.rs;
  halo->np = hf.np;
  memcpy(halo->pos, hf.pos, sizeof(float)*3);
  memcpy(halo->vel, hf.vel, sizeof(float)*3);
  memcpy(halo->J, hf.J, sizeof(float)*3);
  halo->spin = hf.spin;
  halo->phantom = hf.phantom;
  halo->mmp = (hf.flags & MMP_FLAG) ? 1 : 0;
  halo->pid = hf.pid;
  halo->upid = hf.upid;
  halo->orig_id = hf.orig_id;
  halo->last_mm = hf.last_mm;
  halo->desc_pid = -1;
  halo->snapnum = snapnum;
  memcpy(halo->extra_params, hf.extra_params, sizeof(double)*EXTRA_PARAMS);
  return 1;
}

int64_t read_halo_from_file(struct merger_halo *halo, FILE *input, int64_t snapnum) {
  if (!strcasecmp(INPUT_FORMAT, "BINARY"))
    return read_binary_halo_from_file(halo, input, snapnum);

  char buffer[1024];
  while (1) {
    if (!fgets(buffer, 1024, input)) return 0;
    if (buffer[0] != '#') break;
  }
  return read_halo_from_line(halo, buffer, snapnum);
}


int64_t halo_pos_to_tree(struct merger_halo *halo) {
  int64_t idx[3], i;
  for (i=0; i<3; i++) idx[i] = (int64_t)(halo->pos[i]*BOX_DIVISIONS/BOX_WIDTH)
			% (int64_t)(BOX_DIVISIONS);
  return (idx[0]*BOX_DIVISIONS*BOX_DIVISIONS + idx[1]*BOX_DIVISIONS + idx[2]);
}

void *add_to_array(void *array, int64_t *size, int64_t width, void *data) {
  if (!((*size)%1000))
    array = check_realloc(array, width*((*size)+1000), "Allocating array elements.");
  memcpy(array + (*size)*width, data, width);
  *size = (*size) + 1;
  return array;
}

void build_tree(int64_t id, int64_t inputnum) {
  int64_t i = inputnum;
  int64_t id_index = 0;
  double scale, desc_scale;
  struct merger_halo halo;
  int64_t *ids=NULL, num_ids = 0;
  int64_t *new_ids=NULL, num_new_ids = 0;

  ids = add_to_array(ids, &num_ids, sizeof(int64_t), &id_index);

  while (i >= 0) {
    FILE *input = tree_inputs[i];
    scale = output_scales[i];
    desc_scale = (i < num_outputs-1) ? output_scales[i+1] : 0;

    while (1) {
      if (extra_halos[i].id >= 0) {
	halo = extra_halos[i];
	extra_halos[i].id = -1;
      } else {
	if (!read_halo_from_file(&halo, input, i)) break;
      }

      if (halo.descid < 0) {
	fprintf(stderr, "Orphaned halo found!\n");
	exit(1);
      }

      //Check ID of halo to make sure we're still in the right tree:
      for (; id_index < num_ids; id_index++) {
	if (halo.descid == halos[ids[id_index]].id) break;
      }
      if (id_index == num_ids) {
	extra_halos[i] = halo;
	break;
      }
	    
      halo.scale = scale;
      halo.desc_scale = desc_scale;
      new_ids = add_to_array(new_ids, &num_new_ids, sizeof(int64_t), &num_halos);
      halos = add_to_array(halos, &num_halos, sizeof(struct merger_halo), &halo);
    }

    free(ids);
    ids = new_ids;
    num_ids = num_new_ids;
    new_ids = NULL;
    num_new_ids = 0;
    id_index = 0;
    i--;
  }
  free(new_ids);
}


void print_tree_halo(struct merger_halo *h, FILE *output) {
  int64_t next_cop_df = h->next_coprogenitor_df;
  if (next_cop_df > -1) next_cop_df += num_halos_output;
  fprintf(output, " %.4f %8"PRId64" %.4f %8"PRId64" %6"PRId64" %8"PRId64" %8"PRId64" %8"PRId64" %2"PRId64" %.5e %.5e %6f %6f %6f %2"PRId64" %.4f %6f %.5f %.5f %.5f %.3f %.3f %.3f %.3e %.3e %.3e %.5f %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64,
	  h->scale, h->id, h->desc_scale, h->descid, h->num_prog,
	  h->pid, h->upid, h->desc_pid, h->phantom,
	  h->mvir, h->orig_mvir, h->rvir, h->rs, h->vrms,
	  h->mmp, h->last_mm, h->vmax,
	  h->pos[0], h->pos[1], h->pos[2],
	  h->vel[0], h->vel[1], h->vel[2],
	  h->J[0], h->J[1], h->J[2], h->spin,
	  h->breadthfirst_order+num_halos_output, h->depthfirst_order+num_halos_output, h->treeroot_id, h->orig_id,
	  h->snapnum, next_cop_df, h->last_progenitor_df+num_halos_output);
  for (int64_t i=0; i<EXTRA_PARAMS; i++) {
    if (h->extra_params[i] == ((double)((int64_t)h->extra_params[i])))
      fprintf(output, " %"PRId64, (int64_t)h->extra_params[i]);
    else
      fprintf(output, " %g", h->extra_params[i]);
  }
  fprintf(output, "\n");
}

