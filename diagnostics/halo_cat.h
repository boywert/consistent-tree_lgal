#ifndef _HALO_CAT_H_
#define _HALO_CAT_H_

#define CACHE_SIZE 50
#define MASS_LIMIT 3e8

struct halo_data {
  float p[3];
  float v[3];
  float mass;
  float vmax;
  float mass_peak;
  float vmax_acc;
  float vmax_peak;
  float scale;
  float cumul_nd;
  int64_t id, pid, tree_root_id, depthfirst_id, last_prog_id;
  float rvir, vrms;
  struct halo_data *next;
};

struct halo_cat {
  float scale;
  int num_halos, first_halo_index;
  struct halo_data *loc[CACHE_SIZE][CACHE_SIZE][CACHE_SIZE];
  struct halo_cat *next;
  char *filename;
};

extern int box_size;
extern int num_halos;
extern struct halo_data *halos;

void load_halos(char *filename, struct halo_cat **head);
int load_bin_cat(FILE *f);
void load_cat_data(struct halo_cat *cat);
void clear_halos(void);

#endif /* _HALO_CAT_H_ */
