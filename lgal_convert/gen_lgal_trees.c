#include <stdio.h>
#include <inttypes.h>
#include "read_tree.h"
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "gravitational_consistency.h"
#include "gravitational_consistency_vars.h"
#include "halo_io.h"
#include "grav_config.h"
#include "check_syscalls.h"
#include "masses.h"
#include "resort_outputs.h"
#ifndef NO_FORK
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>
#endif /* NO_FORK */
#include "version.h"
#include <math.h>
#include <mpi.h>

struct halo_stash now={0}, evolved = {0};
int64_t *forest_list = NULL;
int64_t num_forests = 0;
float box_size=0;
float max_mvir=0;
float min_mvir=0;
int64_t children = 0;
void do_convert();
int main(int argc, char **argv) {
  int i,j,k,ii,boxd;
  int ierr,rank,size;
  int mark;
  ierr = MPI_Init(&argc, &argv);

  if (argc==1) {
    fprintf(stderr, "Consistent Trees -> LGALAXY Trees, Version %s\n", TREE_VERSION);
    fprintf(stderr, "%s.  See the LICENSE file for redistribution details.\n", TREE_COPYRIGHT);
    fprintf(stderr, "Usage: %s options.cfg\n", argv[0]); exit(1);
  }
  if (argc>1) grav_config(argv[1], 0);

  read_outputs(&(output_scales), &(output_numbers), &(total_outputs));

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
  mark = rank;
  while(mark < BOX_DIVISIONS*BOX_DIVISIONS*BOX_DIVISIONS) {
    i = mark/(BOX_DIVISIONS*BOX_DIVISIONS);
    j = (mark - i*BOX_DIVISIONS*BOX_DIVISIONS)/BOX_DIVISIONS;
    k = (mark - i*BOX_DIVISIONS*BOX_DIVISIONS - j*BOX_DIVISIONS);
    printf("Rank %d: do %d %d %d\n",rank,i,j,k);
    do_convert(i,j,k);
    mark += size;
  }

  ierr = MPI_Finalize();
  return 0;
}
void do_convert(int i, int j, int k) {
  int index = i*BOX_DIVISIONS*BOX_DIVISIONS+j*BOX_DIVISIONS+k;
  char buffer[1024];
  snprintf(buffer,1024,"%s/tree_%d_%d_%d.dat",TREE_OUTBASE,i,j,k);
  printf("%s\n",buffer);
  read_tree(buffer);
  if(all_halos.halos > 0) {
    build_lgal_tree();
    printf("finish building trees\n");
    output_lgal_tree(index);
    printf("%"PRId64" halos found in %s!\n", all_halos.num_halos,buffer);
  }
  else {
    printf("Warning: Ignore %s\n",buffer);
  }
  delete_tree();
}
