#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __APPLE__
#include <mach/mach.h>
#include <mach/mach_time.h>
#else
#include <time.h>
int64_t mach_absolute_time(void) {
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  return t.tv_sec*1e9+t.tv_nsec;
}
#endif /*__APPLE__*/

struct point {
  float pos[3];
  int64_t r, r2, r3, r4;
};

#define FAST3TREE_DIM 3
#define FAST3TREE_TYPE struct point
#include "fast3tree.c"

#define NUM_POINTS 100000000

int main(void) {
  struct point *p = malloc(sizeof(struct point)*NUM_POINTS);
  int64_t i = 0, min_i, max_i, th_num;
  for (i=0; i<NUM_POINTS; i++) {
    p[i].pos[0] = (double)((15*i)%NUM_POINTS)/(double)NUM_POINTS;
    p[i].pos[1] = (double)((122937*i)%NUM_POINTS)/(double)NUM_POINTS;
    p[i].pos[2] = (double)((1939*i)%NUM_POINTS)/(double)NUM_POINTS;
  }
  struct fast3tree *t;
  printf("Starting timer.\n");
  int64_t t1 = mach_absolute_time();
  t = fast3tree_init(NUM_POINTS, p);
  int64_t t2 = mach_absolute_time();
  printf("Time: %.6f s; %p\n", (double)(t2-t1)/1.0e9, t);
  fprintf(stderr, "Verifying...");
  _fast3tree_verify(t->root);
  fprintf(stderr, "ok.\n");
  fast3tree_free(&t);
  t2 = mach_absolute_time();
  t = fast3tree_init(NUM_POINTS, p);
  int64_t t3 = mach_absolute_time();
  printf("Time2: %.6f s; %p\n", (double)(t3-t2)/1.0e9, t);
  fprintf(stderr, "Verifying...");
  _fast3tree_verify(t->root);
  fprintf(stderr, "ok.\n");
  return 0;
}
