#ifndef _ASSEMBLE_HALO_TREES_H_
#define _ASSEMBLE_HALO_TREES_H_
#include <inttypes.h>
#include "merger_halo.h"

void print_tree_halos(int64_t file_id);
void calc_desc_pid_and_num_progs(void);
void calc_last_mm(void);
void calc_mmp(void);
void smooth_mass(void);
int double_sort(const void *a, const void *b);
void conserve_mass(void);
int64_t create_headers(void);
int64_t read_halo_from_file(struct merger_halo *halo, FILE *input, int64_t snapnum);
int64_t halo_pos_to_tree(struct merger_halo *halo);

void *add_to_array(void *array, int64_t *size, int64_t width, void *data);
void build_tree(int64_t id, int64_t inputnum);
void print_tree_halo(struct merger_halo *h, FILE *output);

#endif /* _ASSEMBLE_HALO_TREES_H_ */
