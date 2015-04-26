#include <stdio.h>
#include <inttypes.h>
#include "read_tree.h"

int main(void) {
  read_tree("tree_0_0_0.dat");
  printf("%"PRId64" halos found in tree_0_0_0.dat!\n", all_halos.num_halos);
  return 0;
}
