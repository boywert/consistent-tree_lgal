#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include "check_syscalls.h"
#include "halo_io.h"
#include "tree_halo.h"
#include "grav_config.h"
#include "gravitational_consistency_vars.h"

float box_size, max_mvir, min_mvir;

int64_t is_binary(char *filename) {
  int64_t i;
  FILE *input = check_fopen(filename, "r");
  char buffer[1024];
  int64_t nr = fread(buffer, 1, 1024, input);
  fclose(input);
  if (!nr) return 0;
  for (i=0; i<nr; i++) {
    if (buffer[i] < 1) return 1;
  }
  return 0;
}

void clear_halo_stash(struct halo_stash *h) {
  h->halos = (struct tree_halo *)realloc(h->halos, 0);
  h->id_conv = (int64_t *)realloc(h->id_conv, 0);
  h->max_id = h->min_id = h->num_halos = 0;
}


int main(int argc, char **argv) {
  int64_t i;
  char buffer[1024];
  if (argc < 3) {
    fprintf(stderr, "Usage: %s tree.cfg input1 input2...\n", argv[0]);
    fprintf(stderr, "Converts binary to ascii and vice versa.\n");
  }
  
  grav_config(argv[1], 0);
  
  for (i=2; i<argc; i++) {
    int64_t bin = is_binary(argv[i]);
    snprintf(buffer, 1024, (bin ? "%s.txt" : "%s.bin"), argv[i]);
    INPUT_FORMAT = (bin) ? "BINARY" : "ASCII";
    struct halo_stash h={0};
    fprintf(stderr, (bin) ? "Loading binary %s..." : "Loading ascii %s...", 
	    argv[i]);
    load_halos(argv[i], &h, 0.5, 0);
    fprintf(stderr, "done.\n");
    fprintf(stderr, (bin) ? "Writing ascii %s..." : "Writing binary %s...",
	    buffer);
    INPUT_FORMAT = (bin) ? "ASCII" : "BINARY";
    common_print_halos(buffer, h.halos, h.num_halos, 0);
    clear_halo_stash(&h);
    fprintf(stderr, "done.\n");
  }
  return 0;
}
