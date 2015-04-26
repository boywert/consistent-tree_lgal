#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <memory.h>
#include <math.h>
#include <assert.h>
#include "stringparse.h"
#include "check_syscalls.h"
#include "read_tree.h"

#define SCALE_FACTOR_MUL 10000
struct halo_tree halo_tree = {0};
struct halo_list all_halos = {0};

int sort_by_location(const void *a, const void *b) {
  const struct halo *c = a;
  const struct halo *d = b;
  int64_t dx=c->pos[0]-d->pos[0],
    dy=c->pos[1]-d->pos[1],
    dz=c->pos[2]-d->pos[2];
  if (dx < 0) return -1;
  if (dx > 0) return 1;
  if (dy < 0) return -1;
  if (dy > 0) return 1;
  if (dz < 0) return -1;
  if (dz > 0) return 1;
  return 0;
}

int sort_by_id(const void *a, const void *b) {
  const struct halo_index_key *c = a;
  const struct halo_index_key *d = b;
  if (c->id < d->id) return -1;
  if (d->id < c->id) return 1;
  return 0;
}

int sort_by_desc(const void *a, const void *b) {
  const struct halo *c = a;
  const struct halo *d = b;
  if (c->desc < d->desc) return -1;
  if (d->desc < c->desc) return 1;
  return 0;
}

struct halo *lookup_halo_in_list(struct halo_list *hl, int64_t id) {
  int64_t i;
  int64_t min=0, max=hl->num_halos-1;
  if (id < 0) return 0;
  while (min < max) {
    i = min + ((max-min)/2);
    if (hl->halo_lookup[i].id < id) min = i+1;
    else max = i;
  }
  if (hl->halo_lookup[max].id == id)
    return (&(hl->halos[hl->halo_lookup[max].index]));
  return 0;
}

struct halo_list *lookup_scale(float scale) {
  int64_t index = scale*SCALE_FACTOR_MUL;
  if (index < 0) return 0;
  if (index>=halo_tree.num_scales) return 0;
  index = halo_tree.scale_factor_conv[index];
  if (index < 0) return 0;
  return &(halo_tree.halo_lists[index]);
}

struct halo_list *find_closest_scale(float scale) {
  int64_t index = scale*SCALE_FACTOR_MUL;
  int64_t i=0, j=0;
  if (!halo_tree.num_scales) return NULL;
  if (index < 0) return &(halo_tree.halo_lists[halo_tree.num_lists-1]);
  if (index>=halo_tree.num_scales-1) return halo_tree.halo_lists;
  //index = halo_tree.scale_factor_conv[index];
  while (index+i < halo_tree.num_scales && halo_tree.scale_factor_conv[index+i]<0) i++;
  while (index+j >= 0 && halo_tree.scale_factor_conv[index+j]<0) j--;
  if (index+i < halo_tree.num_scales) {
    index += i;
    if (index-i+j >= 0 && (fabs(j) < fabs(i))) index += j-i;
  } else if (index+j>=0) index += j;
  index = halo_tree.scale_factor_conv[index];
  return &(halo_tree.halo_lists[index]);
}


struct halo_list *_lookup_or_create_scale(float scale) {
  int64_t index = scale*SCALE_FACTOR_MUL;
  int64_t index2;
  assert(index >= 0);
  if (index>=halo_tree.num_scales) {
    int64_t i;
    halo_tree.scale_factor_conv = check_realloc(halo_tree.scale_factor_conv,
						sizeof(int64_t)*(index+1), "Scale factor conversion");
    for (i=halo_tree.num_scales; i<=index; i++) halo_tree.scale_factor_conv[i] = -1;
    halo_tree.num_scales = index+1;
  }
  index2 = halo_tree.scale_factor_conv[index];
  if (index2 < 0) {
    halo_tree.halo_lists = check_realloc(halo_tree.halo_lists,
					 sizeof(struct halo_list)*(halo_tree.num_lists+1),
					 "Halo Lists");
    index2 = halo_tree.num_lists;
    memset(&(halo_tree.halo_lists[index2]), 0, sizeof(struct halo_list));
    halo_tree.halo_lists[index2].scale = scale;
    halo_tree.scale_factor_conv[index] = index2;
    halo_tree.num_lists++;
  }
  return &(halo_tree.halo_lists[index2]);
}


void build_halo_index(struct halo_list *hl) {
  int64_t i;
  hl->halo_lookup = check_realloc(hl->halo_lookup, 
				  sizeof(struct halo_index_key)*hl->num_halos, "Halo lookup index.");
  for (i=0; i<hl->num_halos; i++) {
    hl->halo_lookup[i].id = hl->halos[i].id;
    hl->halo_lookup[i].index = i;
  }
  qsort(hl->halo_lookup, hl->num_halos, sizeof(struct halo_index_key), sort_by_id);
}



void partition_sort_halos(int64_t min, int64_t max,
		struct halo *halos) {
  float minpivot, maxpivot, pivot;
  int64_t i, si;
  struct halo tmp_h;
  if (max-min < 2) return;

  maxpivot = minpivot = halos[min].scale;
  for (i=min+1; i<max; i++) {
    if (halos[i].scale > maxpivot) maxpivot = halos[i].scale;
    if (halos[i].scale < minpivot) minpivot = halos[i].scale;
  }
  if (minpivot==maxpivot) return;
  pivot = minpivot + (maxpivot-minpivot)/2.0;
  si = max-1;
#define SWAP(a,b) {tmp_h = halos[a]; halos[a] = halos[b]; \
    halos[b] = tmp_h; }
  
  for (i=min; i<si; i++)
    if (halos[i].scale <= pivot) { SWAP(i, si); si--; i--; }
  if (i==si && halos[si].scale>pivot) si++;
#undef SWAP
  partition_sort_halos(min, si, halos);
  partition_sort_halos(si, max, halos);
}


void build_tree() {
  int64_t i, j, start, descid;
  struct halo_list *last_hl = 0, *new_hl;
  struct halo *desc;
  memset(&halo_tree, 0, sizeof(struct halo_tree));
  partition_sort_halos(0, all_halos.num_halos, all_halos.halos);
  for (start=0, i=0; i<all_halos.num_halos; i++) {
    if (i>0) assert(all_halos.halos[i].scale <= all_halos.halos[i-1].scale);
    new_hl = _lookup_or_create_scale(all_halos.halos[i].scale);
    if (new_hl != last_hl) {
      new_hl->halos = &(all_halos.halos[i]);
      if (i) {
	last_hl = lookup_scale(all_halos.halos[start].scale);
	last_hl->num_halos = i-start;
      }
      start = i;
      last_hl = new_hl;
    }
  }
  if (last_hl) last_hl->num_halos = i-start;
  if (!halo_tree.num_lists) return;
  last_hl = halo_tree.halo_lists;
  qsort(last_hl->halos, last_hl->num_halos, sizeof(struct halo), sort_by_location);
  build_halo_index(last_hl);
  for (j=0; j<last_hl->num_halos; j++) {
    last_hl->halos[j].parent = lookup_halo_in_list(last_hl, last_hl->halos[j].pid);
    last_hl->halos[j].uparent = lookup_halo_in_list(last_hl, last_hl->halos[j].upid);
    last_hl->halos[j].desc = 0;
  }

  
  for (i=1; i<halo_tree.num_lists; i++) {
    last_hl = &(halo_tree.halo_lists[i-1]);
    new_hl = &(halo_tree.halo_lists[i]);
    for (j=0; j<new_hl->num_halos; j++) {
      descid = (int64_t)(new_hl->halos[j].desc);
      new_hl->halos[j].desc = lookup_halo_in_list(last_hl, (int64_t)descid);
    }
    qsort(new_hl->halos, new_hl->num_halos, sizeof(struct halo), sort_by_desc);
    build_halo_index(new_hl);
    for (j=0; j<new_hl->num_halos; j++) {
      if ((desc = new_hl->halos[j].desc)) {
	if (desc->prog && desc->prog->mvir > new_hl->halos[j].mvir) {
	  new_hl->halos[j].next_coprog = desc->prog->next_coprog;
	  desc->prog->next_coprog = &(new_hl->halos[j]);
	} else {
	  new_hl->halos[j].next_coprog = desc->prog;
	  desc->prog = &(new_hl->halos[j]);
	}
      }
      new_hl->halos[j].parent = lookup_halo_in_list(new_hl, new_hl->halos[j].pid);
      new_hl->halos[j].uparent = lookup_halo_in_list(new_hl, new_hl->halos[j].upid);
    }
  }
}

void read_tree(char *filename) {
  int64_t n;
  FILE *input;
  struct halo h = {0};
  int64_t desc_pid, descid, desc_scale;
  char buffer[1024];

  SHORT_PARSETYPE;
  #define NUM_INPUTS 34
  enum short_parsetype stypes[NUM_INPUTS] = 
    { F, D64, F, D64, D64,    //  #scale id desc_scale desc_id num_prog
      D64, D64, D64, D64,       //   pid upid desc_pid phantom 
      F, F, F, F, F,    //mvir orig_mvir rvir rs vrms 
      D64, F, F,          //mmp? scale_of_last_MM vmax 
      F, F, F, F, F, F,    //x y z vx vy vz
      F, F, F, F, //Jx Jy Jz Spin
      D64, D64, D64, D64, D, D64, D64, //bfid, dfid, trid, ohid, snap, ncdfid, lpdfid
    };
  enum parsetype types[NUM_INPUTS];
  void *data[NUM_INPUTS] = {&(h.scale), &(h.id), &(desc_scale),
                            &(descid), &(h.num_prog), &(h.pid),
			    &(h.upid), &(desc_pid), &(h.phantom), 
			    &(h.mvir), &(h.orig_mvir), &(h.rvir), &(h.rs), &(h.vrms),
			    &(h.mmp), &(h.scale_of_last_MM), &(h.vmax),
			    &(h.pos[0]), &(h.pos[1]), &(h.pos[2]), 
                            &(h.vel[0]), &(h.vel[1]), &(h.vel[2]),
                            &(h.J[0]), &(h.J[1]), &(h.J[2]), &h.spin,
                            &h.breadth_first_id, &h.depth_first_id, &h.tree_root_id, 
                            &h.orig_halo_id, &h.snap_num, &h.next_coprogenitor_depthfirst_id, 
                            &h.last_progenitor_depthfirst_id,
    };

  for (n=0; n<NUM_INPUTS; n++) types[n] = stypes[n];
  input = check_fopen(filename, "r");
  while (fgets(buffer, 1024, input)) {
    if (buffer[0] == '#') continue;
    n = stringparse(buffer, data, (enum parsetype *)types, NUM_INPUTS);
    if (n<NUM_INPUTS) continue;
    h.desc = (struct halo *)(int64_t)descid;
    h.parent = (struct halo *)(int64_t)h.pid;
    h.uparent = (struct halo *)(int64_t)h.upid;
    h.mvir = h.orig_mvir;
    if (!(all_halos.num_halos%3000))
      all_halos.halos = check_realloc(all_halos.halos, sizeof(struct halo)*(all_halos.num_halos+3000), "Allocating Halos.");
   
    all_halos.halos[all_halos.num_halos] = h;
    all_halos.num_halos++;
  }
  fclose(input);

  all_halos.halos = check_realloc(all_halos.halos, sizeof(struct halo)*all_halos.num_halos, "Allocating Halos.");
  build_tree();
}

void delete_tree(void) {
  int64_t i;
  if (!all_halos.num_halos) return;

  for (i=0; i<halo_tree.num_lists; i++)
    free(halo_tree.halo_lists[i].halo_lookup);
  free(halo_tree.halo_lists);
  free(halo_tree.scale_factor_conv);
  memset(&halo_tree, 0, sizeof(struct halo_tree));

  free(all_halos.halos);
  memset(&all_halos, 0, sizeof(struct halo_list));
}

