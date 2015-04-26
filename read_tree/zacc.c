#include <stdio.h>
#include <stdlib.h>
#include "read_tree.h"

float zacc(struct halo *h) {
  if (!h->scale) return 20;
  while (h->prog && (h->pid>=0)) h = h->prog;
  return (1.0/(h->scale) - 1.0);
}

void print_redshift_acc(struct halo *h) {
  float z_a;
  z_a = zacc(h);
  printf("%.5f %d %d %d %d %.3e %.1f %.2f %.2f %.4f %.4f %.4f %.2f %.2f %.2f %f %f\n",
	 h->scale, h->id, h->pid, h->upid, h->phantom, h->mvir,
	 h->rvir, h->rs, h->vmax, h->pos[0], h->pos[1], h->pos[2],
	 h->vel[0], h->vel[1], h->vel[2], h->scale_of_last_MM,
	 z_a);
}


int main(int argc, char **argv) {
  int i, j;
  printf("#Scale ID PID UPID Phantom Mvir Rvir Rs Vmax X Y Z VX VY VZ scale_of_last_mm z_acc\n");
  for (i=1; i<argc; i++) {
    read_tree(argv[i]);
    for (j=0; j<halo_tree.halo_lists[0].num_halos; j++)
      print_redshift_acc(&(halo_tree.halo_lists[0].halos[j]));
    delete_tree();
  }
}
