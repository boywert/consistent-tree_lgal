#include <stdio.h>
#include <stdlib.h>
#include "halo_io.h"
#include "version.h"
#include "gravitational_consistency_vars.h"

struct halo_stash now={0};
float box_size=0;
float max_mvir=0;
float min_mvir=0;

int main(int argc, char **argv) {
  if (argc<2) {
    fprintf(stderr, "Consistent Trees, Version %s\n", TREE_VERSION);
    fprintf(stderr, "(C) 2011-2014, Peter Behroozi.  See the LICENSE file for redistribution details.\n");
    fprintf(stderr, "Usage: %s input.binary\n", argv[0]); exit(1);
  }
  
  char buffer[1024];
  INPUT_FORMAT = "BINARY";
  load_halos(argv[1], &now, 1.0, 0);
  snprintf(buffer, 1024, "%s.ascii", argv[1]);
  INPUT_FORMAT = "ASCII";
  common_print_halos(buffer, now.halos, now.num_halos, 0);
  return 0;
}
