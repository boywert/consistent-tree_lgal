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


struct halo_stash now={0}, evolved = {0};
int64_t *forest_list = NULL;
int64_t num_forests = 0;
float box_size=0;
float max_mvir=0;
float min_mvir=0;
int64_t children = 0;

int main(int argc, char **argv) {
  int i,j,k,findex;
  int64_t total_outputs;
  char buffer[1024];
  if (argc==1) {
    fprintf(stderr, "Consistent Trees -> LGALAXY Trees, Version %s\n", TREE_VERSION);
    fprintf(stderr, "%s.  See the LICENSE file for redistribution details.\n", TREE_COPYRIGHT);
    fprintf(stderr, "Usage: %s options.cfg\n", argv[0]); exit(1);
  }
  if (argc>1) grav_config(argv[1], 0);
  read_outputs(&(output_scales), &(output_numbers), &(total_outputs));
  for(i=0;i<BOX_DIVISIONS;i++)
    for(j=0;j<BOX_DIVISIONS;j++)
      for(k=0;k<BOX_DIVISIONS;k++) {
	findex = i*BOX_DIVISIONS*BOX_DIVISIONS+j*BOX_DIVISIONS+k;
	sprintf(buffer,"tree_%d.dat",(int)i);
	printf("reading %s\n",buffer);
	/* read_tree(treefile); */
	/* build_lgal_tree(); */
	/* output_lgal_tree(findex); */
	/* printf("%"PRId64" halos found in %s!\n", all_halos.num_halos,treefile); */
	/* delete_tree(); */
      }
  return 0;
}
