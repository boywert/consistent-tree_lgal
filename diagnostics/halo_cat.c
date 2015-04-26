#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "distance.h"
#include "halo_cat.h"
#include "stringparse.h"

#define FILENAME_SIZE 250

struct halo_data *halos = NULL;
int num_halos = 0;
int box_size = 0;
float max_p = 0;
float h0 = 0.7;

float readfloat(FILE *input) {
  char buffer[50] = {0};
  float result;
  fgets(buffer, 50, input);
  sscanf(buffer, "%f", &result);
  return result;
}

int sort_vmax_peak(const void *a, const void *b) {
  float vmax_peak_a = halos[*((int *)a)].vmax_peak;
  float vmax_peak_b = halos[*((int *)b)].vmax_peak;
  if (vmax_peak_b < vmax_peak_a) return -1;
  if (vmax_peak_a < vmax_peak_b) return 1;
  return 0;
}

void calc_halo_number_density(struct halo_cat *cat) {
  int *ids = NULL;
  int i;
  float vol = pow(((float)box_size)/h0, 3);
  ids = calloc(cat->num_halos, sizeof(int));
  if (!ids) {
    printf("Couldn't allocate memory to sort ids!\n");
  }
  for (i=0; i<cat->num_halos; i++) ids[i] = cat->first_halo_index+i;
  qsort(ids, cat->num_halos, sizeof(int), sort_vmax_peak);
  for (i=0; i<cat->num_halos; i++)
    halos[ids[i]].cumul_nd = ((float)(i+1))/vol;
  free(ids);
}

void create_halo_caches(struct halo_cat *head) {
  struct halo_cat *cat = head;
  int i, j, max_i, x, y, z;
  for (i=0; i<num_halos; i++) {
    for (j=0; j<3; j++) 
      if (halos[i].p[j] > max_p)
	max_p = halos[i].p[j];
  }
  box_size = max_p + 0.5;

  max_i = cat->first_halo_index + cat->num_halos;
  for (i=cat->first_halo_index; i<max_i; i++) {
    x = halos[i].p[0]*(float)CACHE_SIZE/(float)box_size;
    y = halos[i].p[1]*(float)CACHE_SIZE/(float)box_size;
    z = halos[i].p[2]*(float)CACHE_SIZE/(float)box_size;
    if (x==CACHE_SIZE) x--;
    if (y==CACHE_SIZE) y--;
    if (z==CACHE_SIZE) z--;

    halos[i].next = cat->loc[x][y][z];
    cat->loc[x][y][z] = &(halos[i]);
  }

  calc_halo_number_density(cat);
}

void clear_halos(void) {
  num_halos = 0;
  free(halos);
  halos = NULL;
}

int load_bin_cat(FILE *f) {
  struct stat file_stat;
  int cat_halos = 0;
  fstat(fileno(f), &file_stat);
  if (file_stat.st_size % sizeof(struct halo_data)) return 0;
  
  cat_halos = file_stat.st_size /  sizeof(struct halo_data);
  halos = (struct halo_data *)
    realloc(halos, sizeof(struct halo_data)*(num_halos+cat_halos+1));
  if (!halos) {
    printf("Couldn't malloc room for halos!\n");
    exit(6);
  }

  //Read everything into the halos array.
  fread(&(halos[num_halos]), sizeof(struct halo_data), cat_halos, f);
  num_halos += cat_halos;
  return cat_halos;
}

struct halo_cat *load_cat(char *filename, float scale) {
  struct halo_cat *cat;
  char *dup_filename;

  cat = (struct halo_cat *)malloc(sizeof(struct halo_cat));
  dup_filename = strdup(filename);
  if (!cat || !dup_filename) {
    printf("Couldn't allocate memory for catalog!\n");
    exit(4);
  }
  memset(cat, 0, sizeof(struct halo_cat));

  cat->first_halo_index = num_halos;
  cat->scale = scale;
  cat->filename = dup_filename;
  return cat;
}


void load_cat_data(struct halo_cat *cat) {
  char buffer[1024];
  char bin_fn[FILENAME_SIZE+4];
  int64_t n, descid, np;
  float vmax, rvir, rs, vrms;
  float scale, desc_scale;
  int64_t num_prog, pid, upid, desc_pid, phantom, mmp;
  float orig_mvir, last_mm, mass_acc, vmax_acc;
  FILE *input, *output;
  char *filename;
  struct halo_data h;
  SHORT_PARSETYPE;
  #define NUM_INPUTS 53
  enum short_parsetype stypes[NUM_INPUTS] = 
    { F, D64, F, D64, D64, D64, D64, D64, D64, F, F, F, F, F, D64, F, F, F,  F, F, F,  F, F, K, K, K, K, 
      K, D64, D64, K, K, K, D64, K, K, K, K, K, K, K, K, K, K, K, K, K, K, K,
      F, F, F, F};
  enum parsetype types[NUM_INPUTS];
  void *data[NUM_INPUTS] = 
    {&scale, &(h.id), &desc_scale, &descid, &num_prog, &(h.pid),
     &upid, &desc_pid, &phantom, &(h.mass), &orig_mvir, &(h.rvir), 
     &rs, &(h.vrms), &mmp, &last_mm, &(h.vmax),
     &(h.p[0]), &(h.p[1]), &(h.p[2]), 
     &(h.v[0]), &(h.v[1]), &(h.v[2]),
     NULL, NULL, NULL, NULL, NULL, &(h.depthfirst_id), &(h.tree_root_id),
     NULL, NULL, NULL, &(h.last_prog_id), NULL, NULL, NULL, NULL, NULL, 
     NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
     &mass_acc, &(h.mass_peak), &(h.vmax_acc), &(h.vmax_peak)};


  for (n=0; n<NUM_INPUTS; n++) types[n] = stypes[n];

  filename = cat->filename;
  snprintf(bin_fn, FILENAME_SIZE+4, "%s.bin", filename);
  if ((input = fopen(bin_fn, "rb"))) {
    cat->num_halos = load_bin_cat(input);
    fclose(input);
    create_halo_caches(cat);
    return;
  }

  if (!(input = fopen(filename, "r"))) {
    printf("Could not open halo catalog file %s for reading!\n", filename);
    exit(2);
  }

  h.next = NULL;
  h.scale = cat->scale;
  h.vmax = 0;
  while(fgets(buffer, 1024, input)) {
    n = stringparse(buffer, data, (enum parsetype *)types, NUM_INPUTS);
    if (n!=NUM_INPUTS) continue;

    if (h.mass < MASS_LIMIT) continue;
    if (!(num_halos % 1000)) {
      halos = (struct halo_data *)
	realloc(halos, sizeof(struct halo_data)*(num_halos+1000));
      if (!halos) {
	printf("Couldn't allocate memory for halos!!!\n");
	exit(4);
      }
    }
    halos[num_halos] = h;
    num_halos++;
    cat->num_halos++;
  }
  fclose(input);

  if ((output = fopen(bin_fn, "wb"))) {
    fwrite(&(halos[cat->first_halo_index]), sizeof(struct halo_data), 
	   cat->num_halos, output);
    fclose(output);
  }
  create_halo_caches(cat);
}

void load_halos(char *filename, struct halo_cat **head) {
  FILE *input;
  char buffer[FILENAME_SIZE];
  char halo_file[FILENAME_SIZE];
  float omega_m, omega_l, scale;
  struct halo_cat * new_cat;
  struct halo_cat * last_cat = 0;

  if (!(input = fopen(filename, "rb"))) {
    printf("Could not open halo index file %s for reading!\n", filename);
    exit(2);
  }

  omega_m = readfloat(input);
  omega_l = readfloat(input);
  h0 = readfloat(input);
  if (omega_m <= 0 || omega_m > 1 || omega_l <= 0 || omega_l > 1 ||
      h0 <= 0 || h0 > 1) {
    printf("Incorrect format for halo catalog in %s!\n", filename);
    exit(3);
  }
  if (fabs(omega_m+omega_l-1) > 0.01) {
    printf("Flat cosmologies only.\n");
    exit(4);
  }

  init_cosmology(omega_m, omega_l, h0);
  *head = NULL;
  while (fgets(buffer, FILENAME_SIZE, input)) {
    if (!strlen(buffer)) continue;
    if (sscanf(buffer, "%f %[^\n]", &scale, halo_file) != 2) continue;
    new_cat = load_cat(halo_file, scale);
    if (last_cat) last_cat->next = new_cat;
    else *head = new_cat;
    last_cat = new_cat;
  }
  fclose(input);

}
