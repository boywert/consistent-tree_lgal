#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>

#include "../read_tree/read_tree.h"
#include "../read_tree/read_tree.c"
#include "../read_tree/check_syscalls.h"

void process_tree(float scale,float mass_min, float mass_max);
#define H0 0.7

int main(int argc, char **argv) {
  int64_t i;
  float scale, mass_min, mass_max;
  if (argc<5) {
    printf("Usage: %s scale mass_min mass_max tree.dat ...\n", argv[0]);
    exit(1);
  }

  printf("#Scale Id Mvir X Y Z Scale_Now Id_Now Mvir_now X Y Z\n");
  scale = atof(argv[1]);
  mass_min = atof(argv[2]);
  mass_max = atof(argv[3]);

  for (i=4; i<argc; i++) {
    read_tree(argv[i]);
    process_tree(scale,mass_min,mass_max);
    delete_tree();
  }
  return 0;
}

void process_tree(float scale,float mass_min, float mass_max) {
  struct halo *h;
  struct halo_list *hl = find_closest_scale(scale);
  int64_t i;
  printf("#Scale: %f\n", hl->scale);
  for (i=0; i<hl->num_halos; i++) {
    if (hl->halos[i].mvir < mass_min*H0 ||
	hl->halos[i].mvir > mass_max*H0) continue;
    h = hl->halos + i;
    while (h->desc) h = h->desc;
    printf("%.5f %"PRId64" %.3e %.4f %.4f %.4f %.5f %"PRId64" %.3e %.4f %.4f %.4f\n",
	   hl->halos[i].scale, hl->halos[i].id, hl->halos[i].mvir,
	   hl->halos[i].pos[0], hl->halos[i].pos[1], hl->halos[i].pos[2],
	   h->scale, h->id, h->mvir, h->pos[0], h->pos[1], h->pos[2]);
  }
}

