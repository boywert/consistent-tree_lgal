#ifndef HALO_IO_H
#define HALO_IO_H

#include <stdint.h>
#include "gravitational_consistency.h"

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* From math.h */
#endif /* M_PI */

#define _GC "Gravitational Consistency"
#define _FPAC "Find Parents and Cleanup"
#define _RO "Resort Outputs"
#define _AHT "Assemble Halo Trees"

void gzip_file(char *filename);
void print_halo(FILE *o, struct tree_halo *th);
void common_print_halos(char *filename, struct tree_halo *h, int64_t num_halos, int64_t gzip_halos);
void build_id_conv_list(struct halo_stash *h);
int really_beyond_mmp_ratio(struct tree_halo h1, struct tree_halo h2);
int beyond_mmp_ratio(struct tree_halo h1, struct tree_halo h2);
void add_halo(struct halo_stash *h, struct tree_halo d);
void load_halos(char *filename, struct halo_stash *h, float scale, int dead);
void read_outputs(float **output_scales, int64_t **outputs, int64_t *num_outputs);
void close_timing_log(void);
void print_timing(char *prog, char *info);
#define timed_output(x, ...) { snprintf(buffer, 1024, __VA_ARGS__); print_timing(x, buffer); }

#endif /* HALO_IO_H */
