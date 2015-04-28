#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "grav_config.h"
#include "check_syscalls.h"
#include "resort_outputs.h"
#include <math.h>
#ifndef NO_FORK
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>
#endif /* NO_FORK */
#include "version.h"

#include "stringparse.h"
#include "read_tree.h"
#include "lgal_tree.h"


#define SCALE_FACTOR_MUL 10000
#define GADGET_MASS_CONVERT 1.e-10
struct halo_tree halo_tree = {0};
struct halo_list all_halos = {0};
struct lgal_halo_tree lgal_halo_tree = {0};
struct halo *haloA;
float *output_scales;
int64_t *output_numbers;
int64_t total_outputs;
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

void tree_construct(struct halo *halo, int64_t treenr) {
  struct halo *prog,*next_coprog;
  if(!haloA)
    lgal_halo_tree.root[treenr] = halo;
  else 
    haloA->nexthalo_intree = halo;
  halo->treenr = treenr;
  halo->id_intree = lgal_halo_tree.num_halos_tree[treenr];
  // printf("%" PRId64 " intree_id = %" PRId64 "\n",halo->id,halo->id_intree);
  lgal_halo_tree.num_halos_tree[treenr]++;
  haloA = halo;
  prog = halo->prog;
  if(prog) 
    tree_construct(prog,treenr);
  next_coprog = halo->next_coprog;
  if(next_coprog)
    tree_construct(next_coprog,treenr);
}

long long findlastprogenitorid(struct halo *halo) {
  if(halo->prog)
    if(halo->next_coprog)
      return (long long)halo->next_coprog->id_infile -1;
    else
      return (long long)lgal_halo_tree.num_halos_tree[halo->treenr] - 1;
  else
    return (long long)halo->id_infile;
}

void movetree(int64_t tar, int64_t src) {
  struct halo *cur_halo;
  int64_t count_halo;
  if(tar == src) return;
  count_halo = lgal_halo_tree.num_halos_tree[tar];
  lgal_halo_tree.num_halos_tree[tar] +=  lgal_halo_tree.num_halos_tree[src];
  lgal_halo_tree.num_halos_tree[src] = 0;
  lgal_halo_tree.lastleaf[tar]->nexthalo_intree = lgal_halo_tree.root[src];
  lgal_halo_tree.lastleaf[tar] = lgal_halo_tree.lastleaf[src];
  
  cur_halo = lgal_halo_tree.root[src];
  cur_halo->treenr = tar;
  cur_halo->id_intree = count_halo;
  count_halo++;
  cur_halo = cur_halo->nexthalo_intree;
  while(cur_halo){
    cur_halo->treenr = tar;
    cur_halo->id_intree = count_halo;
    count_halo++;
    cur_halo = cur_halo->nexthalo_intree;
  }
  if(count_halo != lgal_halo_tree.num_halos_tree[tar]) {
    printf("count = %" PRId64 " record = %" PRId64 "\n",count_halo, lgal_halo_tree.num_halos_tree[tar]);
    exit(1);
  }
}

void create_bush(struct halo *halo,int64_t treenr) {
  struct halo *uparent;
  if((uparent = halo->uparent))
    halo = uparent;
  movetree(treenr,halo->treenr);
  halo = halo->nexthalo;
  while(halo) {
    movetree(treenr,halo->treenr);
    halo = halo->nexthalo;
  }
}

void build_lgal_tree() {
  int64_t i,j,count_halo;
  struct halo_list *new_hl;
  struct halo *prog,*cur;
  struct halo *mostmassive,*prev,*temp;
  for (i=2; i<=halo_tree.num_lists; i++) {
    new_hl = &(halo_tree.halo_lists[halo_tree.num_lists-i]);
    for (j=0; j<new_hl->num_halos; j++) {
      if((prog = new_hl->halos[j].prog)) {
	new_hl->halos[j].accu_mass += prog->accu_mass; 
	mostmassive = prog;
      	prog = prog->next_coprog;
      	while(prog) {
	  if(mostmassive->orig_mvir < prog->orig_mvir)
	    mostmassive = prog;
      	  new_hl->halos[j].accu_mass += prog->accu_mass;
	  prog = prog->next_coprog;
      	}
	prog = new_hl->halos[j].prog;
	if(prog != mostmassive) {
	  for(prev = prog; prev != mostmassive; prev = prev->next_coprog);
	  temp = mostmassive->next_coprog;
	  if(prev != prog) 
	    prev->next_coprog = prog;
	  prog = mostmassive;
	  mostmassive->next_coprog = prog->next_coprog;
	  prog->next_coprog = temp;
	}
      }
    }
  }
  new_hl = &(halo_tree.halo_lists[0]);
  lgal_halo_tree.num_trees = new_hl->num_halos;
  lgal_halo_tree.num_halos = all_halos.num_halos;
  lgal_halo_tree.root = check_realloc(lgal_halo_tree.root, new_hl->num_halos*sizeof(struct halo *), "Allocate Tree Root");
  lgal_halo_tree.lastleaf = check_realloc(lgal_halo_tree.lastleaf, new_hl->num_halos*sizeof(struct halo *), "Allocate Tree Last Leaf");
  lgal_halo_tree.num_halos_tree = check_realloc(lgal_halo_tree.num_halos_tree, new_hl->num_halos*sizeof(int64_t), "Allocate Tree Members");
  for(i=0;i<new_hl->num_halos;i++) {
    haloA = 0;
    lgal_halo_tree.num_halos_tree[i] = 0;
    tree_construct(&(new_hl->halos[i]),i);
    lgal_halo_tree.lastleaf[i] = haloA; 
  }
  for(i=0;i<new_hl->num_halos;i++) {
    if(lgal_halo_tree.num_halos_tree[i]) {
        cur = lgal_halo_tree.root[i];
	create_bush(cur, i);
	cur = cur->nexthalo_intree;
	while(cur) {
	  create_bush(cur,i);
	  cur = cur->nexthalo_intree;
	}
    }
  }
  count_halo = 0;
  for(i=0;i<new_hl->num_halos;i++) {
    if(lgal_halo_tree.num_halos_tree[i]) {
      cur = lgal_halo_tree.root[i];
      cur->id_infile = count_halo;
      count_halo++;
      cur = cur->nexthalo_intree;
      while(cur) {
	cur->id_infile = count_halo;
	count_halo++;
	cur = cur->nexthalo_intree;
      }
    }
  }
  if(lgal_halo_tree.num_halos != count_halo) {
    printf("Missing halos\n"); exit(1);
  }
}

int round_to_int( float r ) {
    return floor(r + 0.5); 
}

#define FILENR 0
#define PARTMASS 1.e7

struct lgal_halo_data make_lgal_halo_data(struct halo *halo, int filenr) {
  int i;
  struct lgal_halo_data buffer;
  memset(&buffer,-1,sizeof(struct lgal_halo_data));
  if(halo->desc)
    buffer.Descendant = (int)halo->desc->id_intree;
  if(halo->prog)
    buffer.FirstProgenitor = (int)halo->prog->id_intree;
  if(halo->next_coprog)
    buffer.NextProgenitor = (int)halo->next_coprog->id_intree;
  if(halo->uparent)
    buffer.FirstHaloInFOFgroup = (int)halo->uparent->id_intree;
  else
    buffer.FirstHaloInFOFgroup = (int)halo->id_intree;
  if(halo->nexthalo)
    buffer.NextHaloInFOFgroup = (int)halo->nexthalo->id_intree;

  buffer.Len = (int) round_to_int(halo->orig_mvir/(MASS_RES_OK/1000));
  buffer.M_Mean200 = (float) halo->mvir*GADGET_MASS_CONVERT;
  buffer.M_Crit200 = (float) halo->mvir*GADGET_MASS_CONVERT;
  buffer.M_TopHat = (float) halo->mvir*GADGET_MASS_CONVERT;
  for(i=0;i<3;i++) {
    buffer.Pos[i] = (float) halo->pos[i];
    buffer.Vel[i] = (float) halo->vel[i];
    buffer.Spin[i] = (float) halo->J[i]*GADGET_MASS_CONVERT;
  }
  buffer.VelDisp = (float) halo->vrms;
  buffer.MostBoundID = (long long) 0;
  buffer.SnapNum = (int) output_numbers[halo->snap_num];
  buffer.FileNr = (int) filenr;
  buffer.SubhaloIndex = (int) 0;
  buffer.SubHalfMass = (float) 0;
  return buffer;
}
struct lgal_halo_ids_data make_lgal_halo_ids_data(struct halo *halo, int filenr) {
  struct lgal_halo_ids_data buffer;
  memset(&buffer,-1,sizeof(struct lgal_halo_ids_data));
  buffer.HaloID = (long long) halo->id;
  buffer.FileTreeNr = (long long) filenr;
  if(halo->prog)
    buffer.FirstProgenitor = (long long) halo->prog->id_infile;
  buffer.LastProgenitor = findlastprogenitorid(halo);
  if(halo->nexthalo)
    buffer.NextProgenitor = (long long) halo->nexthalo->id_infile;
  if(halo->desc)
    buffer.Descendant = (long long) halo->desc->id_infile;
  if(halo->uparent)
    buffer.FirstHaloInFOFgroup = (long long) halo->uparent->id_infile;
  else
    buffer.FirstHaloInFOFgroup = (long long) halo->id_infile;
  if(halo->nexthalo)
    buffer.NextHaloInFOFgroup = (long long) halo->nexthalo->id_infile;
#ifdef MAINLEAFID
  buffer.MainLeafID = (long long) 0;
#endif
  buffer.Redshift = (double) (1./output_scales[halo->snap_num] - 1.);
  buffer.PeanoKey = (int)8;
  return buffer;
}
void output_lgal_tree(int filenr) {
  int64_t i,count_tree = 0;
  int buffer;
  char str[1024];
  struct halo *cur;
  struct lgal_halo_data halo_data;
  struct lgal_halo_ids_data halo_ids_data;
  FILE *fp;
  sprintf(str,"mkdir -p %s/treedata",TREE_OUTBASE);
  system(str);
  sprintf(str,"%s/treedata/trees_%03d.%d",TREE_OUTBASE,(int)output_numbers[total_outputs-1],filenr);
  fp = fopen(str,"wb");
  buffer = (int) lgal_halo_tree.num_trees;
  fwrite(&buffer,sizeof(int),1,fp);
  buffer = (int) lgal_halo_tree.num_halos;
  fwrite(&buffer,sizeof(int),1,fp);
  for(i=0;i<lgal_halo_tree.num_trees;i++) {
    if(lgal_halo_tree.num_halos_tree[i]) {
      buffer = lgal_halo_tree.num_halos_tree[i];
      fwrite(&buffer,sizeof(int),1,fp);
      count_tree++;
    }
  }
  for(i=0;i<lgal_halo_tree.num_trees;i++) {
    if(lgal_halo_tree.num_halos_tree[i]) {
      cur = lgal_halo_tree.root[i];
      halo_data = make_lgal_halo_data(cur,filenr);
      fwrite(&halo_data,sizeof(struct lgal_halo_data),1,fp);
      cur = cur->nexthalo_intree;
      while(cur) {
  	halo_data = make_lgal_halo_data(cur,filenr);
  	fwrite(&halo_data,sizeof(struct lgal_halo_data),1,fp);
  	cur = cur->nexthalo_intree;
      }
    }
  }
  rewind(fp);
  lgal_halo_tree.num_trees = count_tree;
  buffer = (int) lgal_halo_tree.num_trees;
  fwrite(&buffer,sizeof(int),1,fp);
  fclose(fp);
  sprintf(str,"%s/treedata/trees_dbids_%03d.%d",TREE_OUTBASE,(int)output_numbers[total_outputs-1],filenr);
  fp = fopen(str,"wb");
  for(i=0;i<lgal_halo_tree.num_trees;i++) {
    if(lgal_halo_tree.num_halos_tree[i]) {
      cur = lgal_halo_tree.root[i];
      halo_ids_data = make_lgal_halo_ids_data(cur,filenr);
      fwrite(&halo_ids_data,sizeof(struct lgal_halo_ids_data),1,fp);
      cur = cur->nexthalo_intree;
      while(cur) {
  	halo_ids_data = make_lgal_halo_ids_data(cur,filenr);
  	fwrite(&halo_ids_data,sizeof(struct lgal_halo_ids_data),1,fp);
  	cur = cur->nexthalo_intree;
      }
    }
  }
  fclose(fp);
}
void build_tree() {
  int64_t i, j, start;
  int check;
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
      new_hl->halos[j].desc = lookup_halo_in_list(last_hl, (int64_t) new_hl->halos[j].descid);
    }
    qsort(new_hl->halos, new_hl->num_halos, sizeof(struct halo), sort_by_desc);
    build_halo_index(new_hl);
    for (j=0; j<new_hl->num_halos; j++) {
      if ((desc = new_hl->halos[j].desc)) {
	new_hl->halos[j].next_coprog = desc->prog;
	desc->prog = &(new_hl->halos[j]);
      }
      new_hl->halos[j].parent = lookup_halo_in_list(new_hl, new_hl->halos[j].pid);
      new_hl->halos[j].uparent = lookup_halo_in_list(new_hl, new_hl->halos[j].upid);
    }
    
    do {
      check = 0;
      for (j=0; j<new_hl->num_halos; j++) {
	if(new_hl->halos[j].uparent) {
	  while(new_hl->halos[j].uparent->uparent){
	    check = 1;
	    printf("1 step up\n");
	    new_hl->halos[j].uparent = new_hl->halos[j].uparent->uparent;
	  }
	}
      }
      for (j=0; j<new_hl->num_halos; j++) {
	if(new_hl->halos[j].uparent) {
	  while(new_hl->halos[j].uparent->uparent){
	    check = 1;
	    printf("2 step up\n");
	    new_hl->halos[j].uparent = new_hl->halos[j].uparent->uparent;
	  }
	}
      }
    } while(check);
    for (j=0; j<new_hl->num_halos; j++) {
      if(new_hl->halos[j].uparent) {
	if(!new_hl->halos[j].uparent->nexthalo)
	  new_hl->halos[j].uparent->nexthalo = &(new_hl->halos[j]);
	else {
	  new_hl->halos[j].nexthalo = new_hl->halos[j].uparent->nexthalo;
	  new_hl->halos[j].uparent->nexthalo = &(new_hl->halos[j]);
	}
      }
    }
  } 
}

void read_tree(char *filename) {
  int64_t n;
  FILE *input;
  struct halo h = {0};
  int64_t desc_pid, desc_scale;
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
                            &(h.descid), &(h.num_prog), &(h.pid),
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
    h.desc = (struct halo *)(int64_t)h.descid;
    h.parent = (struct halo *)(int64_t)h.pid;
    h.uparent = (struct halo *)(int64_t)h.upid;
    h.nexthalo = 0;
    h.nexthalo_intree = 0;
    // h.mvir = h.orig_mvir;
    h.accu_mass = h.mvir;
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

  free(lgal_halo_tree.num_halos_tree);
  free(lgal_halo_tree.root);
  free(lgal_halo_tree.lastleaf);
  memset(&lgal_halo_tree, 0, sizeof(struct lgal_halo_tree));
}

