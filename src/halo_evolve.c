#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include "gravitational_consistency.h"
#include "gravitational_consistency_vars.h"
#include "masses.h"
#include "distance.h"
#include "universe_time.h"
#include "halo_evolve_lib.h"
#include "stringparse.h"
#include "halo_io.h"

struct halo_stash halos = {0};
float box_size = 0; /* Automatically set; in Mpc/h */
float max_mvir = 0; /* Automatically set; in Msun */
float min_mvir = 0;

int main(int argc, char **argv) {
  float a1, a2;
  int64_t i;

  if (argc < 7) {
    printf("Usage: %s scale1 scale2 box_size Om h halolist1 [halolist2 ...]\n", argv[0]);
    exit(1);
  }

  /* Init cosmology and timescales */
  a1 = atof(argv[1]);
  a2 = atof(argv[2]);
  Om = atof(argv[4]);
  if (!Om) Om = 0.27;
  Ol = 1.0 - Om;
  h0 = atof(argv[5]);
  if (!h0) h0 = 0.705;
  init_cosmology(Om, Ol, h0);
  init_time_table(Om, h0);
  gen_ff_cache();

  max_mvir = 0;
  for (i=6; i<argc; i++) load_halos(argv[i], &halos, a1, 0);
  box_size = atof(argv[3]);
  evolve_halos(a1, a2, &halos);
  print_evolved_halos(stdout, &halos);
  return 0;
}
