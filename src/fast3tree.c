/* fast3tree.c
   A fast BSP tree implementation (for arbitrary dimensions).
   (C) 2010-2015, Peter Behroozi, Space Telescope Science Institute.
   Usage:
      #define FAST3TREE_TYPE struct mytype //mytype must have an array pos[]
      #define FAST3TREE_DIM X //Optional---the dimension of pos[], default 3
      #define FAST3TREE_POINTS_PER_LEAF X //Optional, num. points per leaf node
      #define FAST3TREE_PREFIX XYZ //Optional, if using multiple different trees
      #define FAST3TREE_FLOATTYPE float //Optional: or double, or long double
      #include "fast3tree.c"


   PUBLIC METHODS:
   Initialize a fast3tree from a list of points:
      struct fast3tree *fast3tree_init(int64_t n, FAST3TREE_TYPE *p);

   Rebuild a fast3tree from a new (or the same) list of points:
      void fast3tree_rebuild(struct fast3tree *t, int64_t n, FAST3TREE_TYPE *p);

   Rebuilds the tree boundaries, but keeps structure the same:
      void fast3tree_maxmin_rebuild(struct fast3tree *t);

   Frees the tree memory and sets tree pointer to NULL.
      void fast3tree_free(struct fast3tree **t);

   Initialize a fast3tree results structure:
      struct fast3tree_results *fast3tree_results_init(void);

   Find all points within a sphere centered at c[FD] with radius r:
   (FD = FAST3TREE_DIM, usually 3 unless you changed it).
      void fast3tree_find_sphere(struct fast3tree *t,
              struct fast3tree_results *res, float c[FD], float r);

   Find all points within a sphere centered at c[FD] with radius r,
   assuming periodic boundary conditions:
      int fast3tree_find_sphere_periodic(struct fast3tree *t, 
              struct fast3tree_results *res, float c[FD], float r);

   Find all points outside a box with coordinates (x0,y0,...) to (x1,y1,...)
      void fast3tree_find_outside_of_box(struct fast3tree *t, 
              struct fast3tree_results *res, float b[2*FD]);

   Find all points inside a box with coordinates (x0,y0,...) to (x1,y1,...)
      void fast3tree_find_inside_of_box(struct fast3tree *t, 
              struct fast3tree_results *res, float b[2*FD]);

   Reset memory stored in results structure:
      void fast3tree_results_clear(struct fast3tree_results *res);

   Free the results structure returned by fast3tree_find_sphere:
      void fast3tree_results_free(struct fast3tree_results *res);

   END PUBLIC METHODS
*/

#ifndef _FAST3TREE_C_
#define _FAST3TREE_C_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <assert.h>
#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

#ifndef FAST3TREE_PREFIX
#define FAST3TREE_PREFIX
#endif /* FAST3TREE_PREFIX */

#ifndef FAST3TREE_DIM
#define FAST3TREE_DIM 3
#endif /* FAST3TREE_DIM */

#ifndef POINTS_PER_LEAF
#define POINTS_PER_LEAF 40
#endif  /* POINTS_PER_LEAF */

#ifdef FAST3TREE_FLOATTYPE
#define float FAST3TREE_FLOATTYPE
#endif /* FAST3TREE_FLOATTYPE */

#ifndef FAST3TREE_MIN_PARALLEL
#define FAST3TREE_MIN_PARALLEL 2048
#endif /* FAST3TREE_MIN_PARALLEL */

#ifndef FAST3TREE_PARALLEL_CHUNK_SIZE
#define FAST3TREE_PARALLEL_CHUNK_SIZE 200000
#endif /* FAST3TREE_PARALLEL_FACTOR */

#define FAST3TREE_MARKED 1 //Flag for marked nodes

#ifndef FAST3TREE_TYPE
#error Usage:
#error #define FAST3TREE_TYPE point_structure
#error #include "fast3tree.c"
#define FAST3TREE_TYPE struct fast3tree_default_point 
#endif /* FAST3TREE_TYPE */

#define _F3TS(a,b)  a ## b
#define _F3TN(a,b) _F3TS(a, b)
#undef tree3_node
#define tree3_node _F3TN(FAST3TREE_PREFIX,tree3_node)
#undef fast3tree
#define fast3tree  _F3TN(FAST3TREE_PREFIX,fast3tree)
#undef fast3tree_results
#define fast3tree_results _F3TN(FAST3TREE_PREFIX,fast3tree_results)
#undef fast3tree_default_point
#define fast3tree_default_point _F3TN(FAST3TREE_PREFIX,fast3tree_default_point)

struct tree3_node;
struct fast3tree;
struct fast3tree_results;
struct fast3tree_default_point;

#undef fast3tree_init
#define fast3tree_init _F3TN(FAST3TREE_PREFIX,fast3tree_init)
struct fast3tree *fast3tree_init(int64_t n, FAST3TREE_TYPE *p);

#undef fast3tree_rebuild
#define fast3tree_rebuild _F3TN(FAST3TREE_PREFIX,fast3tree_rebuild)
void fast3tree_rebuild(struct fast3tree *t, int64_t n, FAST3TREE_TYPE *p);

#undef fast3tree_maxmin_rebuild
#define fast3tree_maxmin_rebuild _F3TN(FAST3TREE_PREFIX,fast3tree_maxmin_rebuild)
void fast3tree_maxmin_rebuild(struct fast3tree *t);

#undef fast3tree_results_init
#define fast3tree_results_init _F3TN(FAST3TREE_PREFIX,fast3tree_results_init)
struct fast3tree_results *fast3tree_results_init(void);

#undef fast3tree_find_sphere
#define fast3tree_find_sphere _F3TN(FAST3TREE_PREFIX,fast3tree_find_sphere)
static inline void fast3tree_find_sphere(struct fast3tree *t,
			   struct fast3tree_results *res, float c[FAST3TREE_DIM], float r);

#undef fast3tree_find_sphere_skip
#define fast3tree_find_sphere_skip _F3TN(FAST3TREE_PREFIX,fast3tree_find_sphere_skip)
inline void fast3tree_find_sphere_skip(struct fast3tree *t,
		  struct fast3tree_results *res, FAST3TREE_TYPE *tp, float r);

#undef fast3tree_find_sphere_periodic
#define fast3tree_find_sphere_periodic _F3TN(FAST3TREE_PREFIX,fast3tree_find_sphere_periodic)
int fast3tree_find_sphere_periodic(struct fast3tree *t,
			   struct fast3tree_results *res, float c[FAST3TREE_DIM], float r);


#undef fast3tree_find_sphere_marked
#define fast3tree_find_sphere_marked _F3TN(FAST3TREE_PREFIX,fast3tree_find_sphere_marked)
int fast3tree_find_sphere_marked(struct fast3tree *t,
				 struct fast3tree_results *res, float c[FAST3TREE_DIM], float r, int periodic, int do_marking);

#undef fast3tree_find_inside_of_box
#define fast3tree_find_inside_of_box _F3TN(FAST3TREE_PREFIX,fast3tree_find_inside_of_box)
void fast3tree_find_inside_of_box(struct fast3tree *t, struct fast3tree_results *res, float b[2*FAST3TREE_DIM]);

#undef fast3tree_find_outside_of_box
#define fast3tree_find_outside_of_box _F3TN(FAST3TREE_PREFIX,fast3tree_find_outside_of_box)
void fast3tree_find_outside_of_box(struct fast3tree *t, struct fast3tree_results *res, float b[2*FAST3TREE_DIM]);

#undef fast3tree_results_clear
#define fast3tree_results_clear _F3TN(FAST3TREE_PREFIX,fast3tree_results_clear)
void fast3tree_results_clear(struct fast3tree_results *res);

#undef fast3tree_results_free
#define fast3tree_results_free _F3TN(FAST3TREE_PREFIX,fast3tree_results_free)
void fast3tree_results_free(struct fast3tree_results *res);


#ifndef __APPLE__
#ifndef isfinite
#define isfinite(x) finitef(x)
#endif
#endif

#ifdef FAST3TREE_PRINT_TIMING
#define _fast3tree_timing(x) \
  fprintf(stderr, "[Fast3tree Timing] %s: %f s\n", x, (mach_absolute_time()-t0)/1e9); \
  t0 = mach_absolute_time();
#ifdef __APPLE__
#include <mach/mach.h>
#include <mach/mach_time.h>
#else
#include <time.h>
static int64_t mach_absolute_time(void) {
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  return t.tv_sec*1e9+t.tv_nsec;
}
#endif /*__APPLE__*/
#else /* FAST3TREE_PRINT_TIMING */
#define _fast3tree_timing(x)
#endif /* FAST3TREE_PRINT_TIMING */

struct fast3tree_default_point { float pos[FAST3TREE_DIM]; };

struct tree3_node {
  float min[FAST3TREE_DIM], max[FAST3TREE_DIM];
  int64_t num_points;
  int16_t div_dim, flags;
  struct tree3_node *left, *right, *parent;
#ifdef FAST3TREE_EXTRA_INFO
  FAST3TREE_EXTRA_INFO;
#endif /*FAST3TREE_EXTRA_INFO*/
  FAST3TREE_TYPE *points;
};

struct fast3tree {
  FAST3TREE_TYPE *points;
  int64_t num_points;
  struct tree3_node *root;
  int64_t num_nodes;
  int64_t allocated_nodes;
};

struct fast3tree_results {
  int64_t num_points;
  int64_t num_allocated_points;
  FAST3TREE_TYPE **points;
};


/* PRIVATE METHODS */

#undef _fast3tree_check_realloc
#define _fast3tree_check_realloc _F3TN(FAST3TREE_PREFIX,_fast3tree_check_realloc)
void *_fast3tree_check_realloc(void *ptr, size_t size, char *reason);

#define stringify(x) #x
#define to_string(x) stringify(x)
#define _fast3tree_check_realloc_s(x,t,n) x = _fast3tree_check_realloc(x,((int64_t)(n))*((int64_t)(sizeof(t))), "Reallocating " #x " at " __FILE__ ":" to_string(__LINE__));
#define _fast3tree_check_init(x,t,n) t *x = _fast3tree_check_realloc(NULL,((int64_t)(n))*((int64_t)(sizeof(t))), "Allocating " #x " at " __FILE__ ":" to_string(__LINE__));

#undef _fast3tree_build
#define _fast3tree_build _F3TN(FAST3TREE_PREFIX,_fast3tree_build)
void _fast3tree_build(struct fast3tree *t, int64_t within_parallel);

#undef _fast3tree_maxmin_rebuild
#define _fast3tree_maxmin_rebuild _F3TN(FAST3TREE_PREFIX,_fast3tree_maxmin_rebuild)
void _fast3tree_maxmin_rebuild(struct tree3_node *n);



struct fast3tree *fast3tree_init(int64_t n, FAST3TREE_TYPE *p) {
  struct fast3tree *new=NULL;
  new = _fast3tree_check_realloc(new,sizeof(struct fast3tree), "Allocating fast3tree.");
  memset(new, 0, sizeof(struct fast3tree));
  fast3tree_rebuild(new, n, p);
  return new;
}

void fast3tree_rebuild(struct fast3tree *t, int64_t n, FAST3TREE_TYPE *p) {
  t->points = p;
  t->num_points = n;
  _fast3tree_build(t, 0);
}

void fast3tree_maxmin_rebuild(struct fast3tree *t) {
  _fast3tree_maxmin_rebuild(t->root);
}

#undef fast3tree_free
#define fast3tree_free _F3TN(FAST3TREE_PREFIX,fast3tree_free)
void fast3tree_free(struct fast3tree **t) {
  if (!t) return;
  struct fast3tree *u = *t;
  if (u) {
    free(u->root);
    free(u);
  }
  *t = NULL;
}

#undef _fast3tree_box_inside_box
#define _fast3tree_box_inside_box \
  _F3TN(FAST3TREE_PREFIX,_fast3tree_box_inside_box)
static inline int _fast3tree_box_inside_box(struct tree3_node *node, float b[2*FAST3TREE_DIM]) {
  int i;
  for (i=0; i<FAST3TREE_DIM; i++) {
    if (node->max[i]>b[i+FAST3TREE_DIM]) return 0;
    if (node->min[i]<b[i]) return 0;
  }
  return 1;
}


#undef _fast3tree_box_intersect_box
#define _fast3tree_box_intersect_box \
  _F3TN(FAST3TREE_PREFIX,_fast3tree_box_intersect_box)
static inline int _fast3tree_box_intersect_box(struct tree3_node *node, float b[2*FAST3TREE_DIM]) {
  int i;
  for (i=0; i<FAST3TREE_DIM; i++) {
    if (node->max[i]<b[i]) return 0;
    if (node->min[i]>b[i+FAST3TREE_DIM]) return 0;
  }
  return 1;
}


#undef _fast3tree_box_not_intersect_sphere
#define _fast3tree_box_not_intersect_sphere \
  _F3TN(FAST3TREE_PREFIX,_fast3tree_box_not_intersect_sphere)
/* Accurate, fast */
static inline int _fast3tree_box_not_intersect_sphere(const struct tree3_node *node, const float c[FAST3TREE_DIM], const float r) {
  int i;
  float d = 0, e;
  const float r2 = r*r;
  for (i=0; i<FAST3TREE_DIM; i++) {
    if ((e = c[i]-node->min[i])<0) {
      d+=e*e;
      if (d >= r2) return 1;
    } else if ((e = c[i]-node->max[i])>0) {
      d+=e*e;
      if (d >= r2) return 1;
    }
  }
  return 0;
}

/* Fast, accurate. */
#undef _fast3tree_box_inside_sphere
#define _fast3tree_box_inside_sphere _F3TN(FAST3TREE_PREFIX,_fast3tree_box_inside_sphere)
inline int _fast3tree_box_inside_sphere(const struct tree3_node *node, const float c[FAST3TREE_DIM], const float r) {
  int i;
  float dx, dx2, dist = 0, r2 = r*r;
  if (fabs(c[0]-node->min[0]) > r) return 0; //Rapid short-circuit.
  for (i=0; i<FAST3TREE_DIM; i++) {
    dx = node->min[i] - c[i];
    dx *= dx;
    dx2 = c[i]-node->max[i];
    dx2 *= dx2;
    if (dx2 > dx) dx = dx2;
    dist += dx;
    if (dist > r2) return 0;
  }
  return 1;
}

#undef _fast3tree_sphere_inside_box
#define _fast3tree_sphere_inside_box _F3TN(FAST3TREE_PREFIX,_fast3tree_sphere_inside_box)
static inline int _fast3tree_sphere_inside_box(const struct tree3_node *node, const float c[FAST3TREE_DIM], const float r) {
  int i;
  for (i=0; i<FAST3TREE_DIM; i++) {
    if (c[i]-r < node->min[i]) return 0;
    if (c[i]+r > node->max[i]) return 0;
  }
  return 1;
}

#undef _fast3tree_check_results_space
#define _fast3tree_check_results_space _F3TN(FAST3TREE_PREFIX,_fast3tree_check_results_space)
static inline void _fast3tree_check_results_space(const struct tree3_node *n, struct fast3tree_results *res) {
  if (res->num_points + n->num_points > res->num_allocated_points) {
    res->num_allocated_points = res->num_points + n->num_points + 1000;
    res->points = _fast3tree_check_realloc(res->points, 
 res->num_allocated_points * sizeof(FAST3TREE_TYPE *), "Allocating fast3tree results");
  }
}


#undef _fast3tree_mark_node
#define _fast3tree_mark_node _F3TN(FAST3TREE_PREFIX,_fast3tree_mark_node)
static inline void _fast3tree_mark_node(struct tree3_node *n, int16_t mark) {
  n->flags = mark;
  if (n->div_dim > -1) {
    _fast3tree_mark_node(n->left, mark);
    _fast3tree_mark_node(n->right, mark);
  }
}

#undef _fast3tree_find_sphere
#define _fast3tree_find_sphere _F3TN(FAST3TREE_PREFIX,_fast3tree_find_sphere)
void _fast3tree_find_sphere(const struct tree3_node *n, struct fast3tree_results *res, float c[FAST3TREE_DIM], const float r) {
  int64_t i,j;
  float r2, dist, dx;

  if (_fast3tree_box_not_intersect_sphere(n,c,r)) return;
#if FAST3TREE_DIM < 6
  if (_fast3tree_box_inside_sphere(n,c,r)) { /* Entirely inside sphere */
    _fast3tree_check_results_space(n,res);
    for (i=0; i<n->num_points; i++)
      res->points[res->num_points+i] = n->points+i;
    res->num_points += n->num_points;
    return;
  }
#endif /* FAST3TREE_DIM < 6 */

  if (n->div_dim < 0) { /* Leaf node */
    r2 = r*r;
    _fast3tree_check_results_space(n,res);
    for (i=0; i<n->num_points; i++) {
      j = dist = 0;
      float *pos = n->points[i].pos;
      for (; j<FAST3TREE_DIM; j++) {
	dx = c[j]-pos[j];
	dist += dx*dx;
      }
      if (dist < r2) {
	res->points[res->num_points] = n->points + i;
	res->num_points++;
      }
    }
    return;
  }
  _fast3tree_find_sphere(n->left, res, c, r);
  _fast3tree_find_sphere(n->right, res, c, r);
}


#undef _fast3tree_find_sphere_skip
#define _fast3tree_find_sphere_skip _F3TN(FAST3TREE_PREFIX,_fast3tree_find_sphere_skip)
void _fast3tree_find_sphere_skip(const struct tree3_node *n, struct fast3tree_results *res, float c[FAST3TREE_DIM], const float r, FAST3TREE_TYPE *tp) {
  int64_t i,j;
  float r2, dist, dx;

  if (n->points + n->num_points <= tp) return; //Skip already-processed nodes
  if (_fast3tree_box_not_intersect_sphere(n,c,r)) return;
#if FAST3TREE_DIM < 6
  if (_fast3tree_box_inside_sphere(n,c,r)) { /* Entirely inside sphere */
    _fast3tree_check_results_space(n,res);
    for (i=0; i<n->num_points; i++)
      res->points[res->num_points+i] = n->points+i;
    res->num_points += n->num_points;
    return;
  }
#endif /* FAST3TREE_DIM < 6 */

  if (n->div_dim < 0) { /* Leaf node */
    r2 = r*r;
    _fast3tree_check_results_space(n,res);
    i=0;
    if (n->points < tp) {
      res->points[res->num_points] = tp;
      res->num_points++;
      i = (tp-n->points)+1;
    }
    for (; i<n->num_points; i++) {
      j = dist = 0;
      float *pos = n->points[i].pos;
      for (; j<FAST3TREE_DIM; j++) {
	dx = c[j]-pos[j];
	dist += dx*dx;
      }
      if (dist < r2) {
	res->points[res->num_points] = n->points + i;
	res->num_points++;
      }
    }
    return;
  }
  _fast3tree_find_sphere_skip(n->left, res, c, r, tp);
  _fast3tree_find_sphere_skip(n->right, res, c, r, tp);
}


struct fast3tree_results *fast3tree_results_init(void) {
  struct fast3tree_results *res = 
    _fast3tree_check_realloc(NULL, sizeof(struct fast3tree_results), "Allocating fast3tree results structure.");
  res->points = NULL;
  res->num_points = 0;
  res->num_allocated_points = 0;
  return res;
}

static inline void fast3tree_find_sphere(struct fast3tree *t, struct fast3tree_results *res, float c[FAST3TREE_DIM], float r) {
  res->num_points = 0;
  _fast3tree_find_sphere(t->root, res, c, r);
}

inline void fast3tree_find_sphere_skip(struct fast3tree *t, struct fast3tree_results *res, FAST3TREE_TYPE *tp, float r) {
  res->num_points = 0;
  _fast3tree_find_sphere_skip(t->root, res, tp->pos, r, tp);
}


/* Guaranteed to be stable with respect to floating point round-off errors.*/
#undef _fast3tree_find_sphere_offset
#define _fast3tree_find_sphere_offset _F3TN(FAST3TREE_PREFIX,_fast3tree_find_sphere_offset)
void _fast3tree_find_sphere_offset(struct tree3_node *n, struct fast3tree_results *res, float c[FAST3TREE_DIM], float c2[FAST3TREE_DIM], float o[FAST3TREE_DIM], const float r, const int marked, const int do_marking) {
  int64_t i,j;
  float r2, dist, dx;

  int64_t onlyone = (marked && (n->flags & FAST3TREE_MARKED)) ? 1 : 0;

  if (_fast3tree_box_not_intersect_sphere(n,c2,r*1.01)) return;
  if (_fast3tree_box_inside_sphere(n,c2,r*0.99)) { /* Entirely inside sphere */
    if (do_marking) n->flags |= FAST3TREE_MARKED;
    _fast3tree_check_results_space(n,res);
    if (onlyone) {
      res->points[res->num_points++] = n->points;
      return;
    }
    for (i=0; i<n->num_points; i++)
      res->points[res->num_points+i] = n->points+i;
    res->num_points += n->num_points;
    return;
  }

  if (n->div_dim < 0) { /* Leaf node */
    r2 = r*r;
    _fast3tree_check_results_space(n,res);
    for (i=0; i<n->num_points; i++) {
      j = dist = 0;
      float *pos = n->points[i].pos;
      for (; j<FAST3TREE_DIM; j++) {
	dx = o[j] - fabs(c[j]-pos[j]);
	dist += dx*dx;
      }
      if (dist < r2) {
	res->points[res->num_points] = n->points + i;
	res->num_points++;
	if (onlyone) return;
      }
    }
    return;
  }
  int64_t cur_points = res->num_points;
  _fast3tree_find_sphere_offset(n->left, res, c, c2, o, r, marked, do_marking);
  if (onlyone && (cur_points < res->num_points)) return;
  _fast3tree_find_sphere_offset(n->right, res, c, c2, o, r, marked, do_marking);
}

#undef _fast3tree_find_sphere_periodic_dim
#define _fast3tree_find_sphere_periodic_dim _F3TN(FAST3TREE_PREFIX,_fast3tree_find_sphere_periodic_dim)
void _fast3tree_find_sphere_periodic_dim(struct fast3tree *t, struct fast3tree_results *res, float c[FAST3TREE_DIM], float c2[FAST3TREE_DIM], float o[FAST3TREE_DIM], float r, float dims[FAST3TREE_DIM], int dim, int marked, int do_marking) {
  float c3[FAST3TREE_DIM];
  if (dim<0) {
    _fast3tree_find_sphere_offset(t->root, res, c, c2, o, r, marked, do_marking);
    return;
  }
  memcpy(c3, c2, sizeof(float)*FAST3TREE_DIM);
  o[dim]=0;
  _fast3tree_find_sphere_periodic_dim(t, res, c, c3, o, r, dims, dim-1, marked, do_marking);
  if (c[dim]+r > t->root->max[dim]) {
    c3[dim] = c[dim]-dims[dim];
    o[dim] = dims[dim];
    _fast3tree_find_sphere_periodic_dim(t, res, c, c3, o, r, dims, dim-1, marked, do_marking);
  }
  if (c[dim]-r < t->root->min[dim]) {
    c3[dim] = c[dim]+dims[dim];
    o[dim] = dims[dim];
    _fast3tree_find_sphere_periodic_dim(t, res, c, c3, o, r, dims, dim-1, marked, do_marking);
  }
}

int fast3tree_find_sphere_periodic(struct fast3tree *t, struct fast3tree_results *res, float c[FAST3TREE_DIM], float r) {
  float dims[FAST3TREE_DIM], o[FAST3TREE_DIM];
  int i;
  
  if (_fast3tree_sphere_inside_box(t->root, c, r)) {
    fast3tree_find_sphere(t, res, c, r);
    return 2;
  }

  for (i=0; i<FAST3TREE_DIM; i++) {
    dims[i] = t->root->max[i] - t->root->min[i];
    o[i] = 0;
    if (r*2.0 > dims[i]) return 0; //Avoid wraparound intersections.
  }

  res->num_points = 0;
  _fast3tree_find_sphere_periodic_dim(t, res, c, c, o, r, dims, FAST3TREE_DIM-1, 0, 0);
  return 1;
}


int fast3tree_find_sphere_marked(struct fast3tree *t, struct fast3tree_results *res, float c[FAST3TREE_DIM], float r, int periodic, int do_marking) {
  float dims[FAST3TREE_DIM], o[FAST3TREE_DIM] = {0};
  int i;
  
  res->num_points = 0;
  if (!periodic || _fast3tree_sphere_inside_box(t->root, c, r)) {
     _fast3tree_find_sphere_offset(t->root, res, c, c, o, r, 1, do_marking);
    return 2;
  }

  for (i=0; i<FAST3TREE_DIM; i++) {
    dims[i] = t->root->max[i] - t->root->min[i];
    o[i] = 0;
    if (r*2.0 > dims[i]) return 0; //Avoid wraparound intersections.
  }

  _fast3tree_find_sphere_periodic_dim(t, res, c, c, o, r, dims, FAST3TREE_DIM-1, 1, do_marking);
  return 1;
}

#undef _fast3tree_find_outside_of_box
#define _fast3tree_find_outside_of_box _F3TN(FAST3TREE_PREFIX,_fast3tree_find_outside_of_box)
static inline void _fast3tree_find_outside_of_box(struct tree3_node *n, struct fast3tree_results *res, float b[2*FAST3TREE_DIM]) {
  int64_t i,j;
  if (_fast3tree_box_inside_box(n, b)) return;
  if (!_fast3tree_box_intersect_box(n, b)) {
    _fast3tree_check_results_space(n,res);
    for (i=0; i<n->num_points; i++) {
      res->points[res->num_points] = n->points+i;
      res->num_points++;
    }
  }
  else {
    if (n->div_dim < 0) {
      _fast3tree_check_results_space(n,res);
      for (i=0; i<n->num_points; i++) {
	for (j=0; j<FAST3TREE_DIM; j++) {
	  if (n->points[i].pos[j] < b[j]) break;
	  if (n->points[i].pos[j] > b[j+FAST3TREE_DIM]) break;
	}
	if (j==FAST3TREE_DIM) continue; //Inside
	res->points[res->num_points] = n->points+i;
	res->num_points++;
      }
    }
    else {
      _fast3tree_find_outside_of_box(n->left, res, b);
      _fast3tree_find_outside_of_box(n->right, res, b);
    }
  }
}

void fast3tree_find_outside_of_box(struct fast3tree *t, struct fast3tree_results *res, float b[2*FAST3TREE_DIM]) {
  res->num_points = 0;
  _fast3tree_find_outside_of_box(t->root, res, b);
}

#undef _fast3tree_find_inside_of_box
#define _fast3tree_find_inside_of_box _F3TN(FAST3TREE_PREFIX,_fast3tree_find_inside_of_box)
static inline void _fast3tree_find_inside_of_box(struct tree3_node *n, struct fast3tree_results *res, float b[2*FAST3TREE_DIM]) {
  int64_t i,j;
  if (!_fast3tree_box_intersect_box(n, b)) return;
  if (_fast3tree_box_inside_box(n, b)) {
    _fast3tree_check_results_space(n,res);
    for (i=0; i<n->num_points; i++) {
      res->points[res->num_points] = n->points+i;
      res->num_points++;
    }
  }
  else {
    if (n->div_dim < 0) {
      _fast3tree_check_results_space(n,res);
      for (i=0; i<n->num_points; i++) {
	for (j=0; j<FAST3TREE_DIM; j++) {
	  if (n->points[i].pos[j] < b[j]) break;
	  if (n->points[i].pos[j] > b[j+FAST3TREE_DIM]) break;
	}
	if (j!=FAST3TREE_DIM) continue;
	res->points[res->num_points] = n->points+i;
	res->num_points++;
      }
    }
    else {
      _fast3tree_find_inside_of_box(n->left, res, b);
      _fast3tree_find_inside_of_box(n->right, res, b);
    }
  }
}

void fast3tree_find_inside_of_box(struct fast3tree *t, struct fast3tree_results *res, float b[2*FAST3TREE_DIM]) {
  res->num_points = 0;
  _fast3tree_find_inside_of_box(t->root, res, b);
}


void fast3tree_results_clear(struct fast3tree_results *res) {
  if (res->points) free(res->points);
  memset(res, 0, sizeof(struct fast3tree_results));
}


void fast3tree_results_free(struct fast3tree_results *res) {
  if (!res) return;
  if (res->points) free(res->points);
  memset(res, 0, sizeof(struct fast3tree_results));
  free(res);
}


#undef _fast3tree_find_largest_dim
#define _fast3tree_find_largest_dim _F3TN(FAST3TREE_PREFIX,_fast3tree_find_largest_dim)
static inline int64_t _fast3tree_find_largest_dim(float * min, float * max) {
  int64_t i, dim = FAST3TREE_DIM-1;
  float d = max[FAST3TREE_DIM-1]-min[FAST3TREE_DIM-1], d2;
  for (i=0; i<(FAST3TREE_DIM-1); i++) {
    d2 = max[i] - min[i];
    if (d2 > d) { d=d2; dim = i; }
  }
  return dim;
}

#undef _fast3tree_sort_dim_pos
#define _fast3tree_sort_dim_pos _F3TN(FAST3TREE_PREFIX,_fast3tree_sort_dim_pos)
static inline int64_t _fast3tree_sort_dim_pos(struct tree3_node * node,
					      float balance) {
  int64_t dim = node->div_dim = 
    _fast3tree_find_largest_dim(node->min, node->max);
  FAST3TREE_TYPE *p = node->points;
  int64_t i,j=node->num_points-1;
  FAST3TREE_TYPE temp;
  float lim = 0.5*(node->max[dim]+node->min[dim]);

  if (node->max[dim]==node->min[dim]) return node->num_points;

  for (i=0; i<j; i++) {
    if (p[i].pos[dim] > lim) {
      temp = p[j];
      p[j] = p[i];
      p[i] = temp;
      i--;
      j--;
    }
  }
  if ((i==j) && (p[i].pos[dim] <= lim)) i++;
  return i;
}

#undef _fast3tree_find_minmax
#define _fast3tree_find_minmax _F3TN(FAST3TREE_PREFIX,_fast3tree_find_minmax)
static inline void _fast3tree_find_minmax(struct tree3_node *node) {
  int64_t i, j;
  float x;
  FAST3TREE_TYPE * p = node->points;
  assert(node->num_points > 0);
  for (j=0; j<FAST3TREE_DIM; j++) node->min[j] = node->max[j] = p[0].pos[j];
  for (i=1; i<node->num_points; i++)  {
    for (j=0; j<FAST3TREE_DIM; j++) {
      x = p[i].pos[j];
      if (x<node->min[j]) node->min[j] = x;
      else if (x>node->max[j]) node->max[j] = x;
    }
  }
}

#undef _fast3tree_verify
#define _fast3tree_verify _F3TN(FAST3TREE_PREFIX,_fast3tree_verify)
void _fast3tree_verify(struct tree3_node *node) {
  struct tree3_node n = *node;
  int64_t i;
  if (node->num_points>0)
    _fast3tree_find_minmax(&n);
  for (i=0; i<FAST3TREE_DIM; i++) {
    if (n.min[i] < node->min[i]) {
      struct tree3_node *parent = node->parent;
      while (parent != parent->parent) parent = parent->parent;
      fprintf(stderr, "%ld: exp min: %.7f; actual min: %.7f\n", (node-parent), node->min[i], n.min[i]);
    }
    assert(n.min[i] >= node->min[i]);
    assert(n.max[i] <= node->max[i]);
  }
  if (n.div_dim > -1) {
    assert(node->left->max[n.div_dim] < node->right->max[n.div_dim]);
    _fast3tree_verify(node->left);
    _fast3tree_verify(node->right);
  }
}

#undef _fast3tree_split_node
#define _fast3tree_split_node _F3TN(FAST3TREE_PREFIX,_fast3tree_split_node)
void _fast3tree_split_node(struct fast3tree *t, struct tree3_node *node) {
  int64_t num_left;
  struct tree3_node *null_ptr = NULL;
  struct tree3_node *left, *right;
  int64_t left_index, node_index;

  num_left = _fast3tree_sort_dim_pos(node, 1.0);
  if (num_left == node->num_points || num_left == 0) 
  { //In case all node points are at same spot
    node->div_dim = -1;
    return;
  }

  node_index = node - t->root;
  if ((t->num_nodes+2) > t->allocated_nodes) {
    t->allocated_nodes = t->allocated_nodes*1.05 + 1000;
    t->root = _fast3tree_check_realloc(t->root, sizeof(struct tree3_node)*(t->allocated_nodes), "Tree nodes");
    node = t->root + node_index;
  }

  left_index = t->num_nodes;
  t->num_nodes+=2;

  node->left = null_ptr + left_index;
  node->right = null_ptr + (left_index + 1);

  left = t->root + left_index;
  right = t->root + (left_index + 1);
  memset(left, 0, sizeof(struct tree3_node)*2);

  right->parent = left->parent = null_ptr + node_index;
  left->num_points = num_left;
  right->num_points = node->num_points - num_left;
  left->div_dim = right->div_dim = -1;
  left->points = node->points;
  right->points = node->points + num_left;

  _fast3tree_find_minmax(left);
  _fast3tree_find_minmax(right);

  if (left->num_points > POINTS_PER_LEAF)
    _fast3tree_split_node(t, left);

  right = t->root + (left_index + 1);
  if (right->num_points > POINTS_PER_LEAF)
    _fast3tree_split_node(t, right);
}

#ifdef _OPENMP
#undef _fast3tree_preprocess_parallel
#define _fast3tree_preprocess_parallel _F3TN(FAST3TREE_PREFIX,_fast3tree_preprocess_parallel)
static inline void _fast3tree_preprocess_parallel(struct fast3tree *t, struct tree3_node *node) {
  int64_t i, j, k;
  FAST3TREE_TYPE * p = node->points;
  assert(node->num_points > 0);
  int64_t cpus;
#pragma omp parallel
  cpus = omp_get_num_threads();

  _fast3tree_check_init(min_i, int64_t, cpus);
  _fast3tree_check_init(max_i, int64_t, cpus);
  _fast3tree_check_init(min, float, cpus*FAST3TREE_DIM);
  _fast3tree_check_init(max, float, cpus*FAST3TREE_DIM);

#ifdef FAST3TREE_PRINT_TIMING
  int64_t t0 = mach_absolute_time();
#endif /* FAST3TREE_PRINT_TIMING */
  int64_t block_size = ceil((float)node->num_points/(float)cpus);

#pragma omp parallel default(none) shared(max_i,min_i,min,max,block_size,node,p)
  {
    int64_t tid = omp_get_thread_num();
    int64_t toff = tid*FAST3TREE_DIM;
    int64_t mii = tid*block_size;
    int64_t mai = mii + block_size;
    float my_min[FAST3TREE_DIM]={0}, my_max[FAST3TREE_DIM]={0};
    int64_t i,j;
    if (mai > node->num_points) mai = node->num_points;
    for (i=0; i<node->num_points; i++) {
      for (j=0; j<FAST3TREE_DIM; j++) {
	if (!isfinite(p[i].pos[j])) break;
	my_min[j] = my_max[j] = p[i].pos[j];
      }
      if (j==FAST3TREE_DIM) break;
    }
    for (i=mii; i<mai; i++) {
      for (j=0; j<FAST3TREE_DIM; j++)
	if (!isfinite(p[i].pos[j])) break;
      if (j<FAST3TREE_DIM) {
	FAST3TREE_TYPE tmp = p[i];
	mai--;
	p[i] = p[mai];
	p[mai] = tmp;
	i--;
	continue;
      }
      for (j=0; j<FAST3TREE_DIM; j++) {
	float x = p[i].pos[j];
	if (x<my_min[j]) my_min[j] = x;
	else if (x>my_max[j]) my_max[j] = x;
      }
    }

    max_i[tid] = mai;
    min_i[tid] = mii;
    memcpy(max+toff, my_max, sizeof(float)*FAST3TREE_DIM);
    memcpy(min+toff, my_min, sizeof(float)*FAST3TREE_DIM);
  }

  _fast3tree_timing("Parallel min/max+finite-check time");

  for (j=0; j<FAST3TREE_DIM; j++) {
    node->max[j] = max[j];
    node->min[j] = min[j];
    for (i=1; i<cpus; i++) {
      if (node->max[j] < max[j+i*FAST3TREE_DIM])
	node->max[j] = max[j+i*FAST3TREE_DIM];
      if (node->min[j] > min[j+i*FAST3TREE_DIM])
	node->min[j] = min[j+i*FAST3TREE_DIM];
    }
  }

  //Reorder any infinite points
  j = cpus-1;
  for (i=0; i<j; i++) {
    int64_t exp_max = min_i[i]+block_size;
    int64_t to_move = exp_max - max_i[i];
    while (j>i) {
      int64_t avail = max_i[j]-min_i[j];
      int64_t to_swap = avail;
      if (to_swap > to_move) to_swap = to_move;
      for (k = 0; k<to_swap; k++) {
	max_i[j]--;
	FAST3TREE_TYPE tmp = p[max_i[i]];
	p[max_i[i]] = p[max_i[j]];
	p[max_i[j]] = tmp;
	max_i[i]++;
      }
      if (to_swap == avail) j--;
      else break;
    }
  }

  //Recalculate number of points
  node->num_points = 0;
  for (i=0; i<cpus; i++) node->num_points += max_i[i]-min_i[i];
  t->num_points = node->num_points;
  if (!node->num_points) {
    memset(node->min, 0, sizeof(float)*FAST3TREE_DIM);
    memset(node->max, 0, sizeof(float)*FAST3TREE_DIM);
  }

  free(max);
  free(min);
  free(max_i);
  free(min_i);
  _fast3tree_timing("Cleanup");
}


#undef _fast3tree_parallel_mapping
#define _fast3tree_parallel_mapping _F3TN(FAST3TREE_PREFIX,_fast3tree_parallel_mapping)
int64_t _fast3tree_parallel_mapping(struct fast3tree *t, struct tree3_node *node, int64_t *factors, int64_t *factors_prod, float *chunk_size, int64_t *mapping, int64_t last_div_dim, int64_t last_block) {
  int64_t i, n, div_dim, prev_factor;
  struct tree3_node *null_ptr = NULL;
  for (i=0; i<FAST3TREE_DIM; i++) {
    div_dim = (i+1+last_div_dim) % FAST3TREE_DIM;
    if (factors[div_dim] > 1) break;
  }
  if (factors[div_dim] > 1) {
    prev_factor = factors[div_dim];
    int64_t left_chunks = factors[div_dim]/2;
    int64_t right_chunks = factors[div_dim] - left_chunks;
    int64_t node_index = node - t->root;
    if ((t->num_nodes+2) > t->allocated_nodes) {
      t->allocated_nodes = t->allocated_nodes*1.05 + 1000;
      t->root = _fast3tree_check_realloc(t->root, sizeof(struct tree3_node)*(t->allocated_nodes), "Tree nodes");
      node = t->root + node_index;
    }
    
    node->div_dim = div_dim;
    int64_t left_index = t->num_nodes;
    t->num_nodes+=2;

    node->left = null_ptr + left_index;
    node->right = null_ptr + (left_index + 1);

    struct tree3_node *left = t->root + left_index;
    struct tree3_node *right = t->root + (left_index + 1);
    *left = *node;
    *right = *node;
    left->left = left->right = right->left = right->right = 0;
    left->div_dim = right->div_dim = -1;
    right->parent = left->parent = null_ptr + node_index;

    //Need to calculate boundaries stably wrt. floating-point errors
    int64_t total_chunks = left_chunks + ((float)(left->min[div_dim]-t->root->min[div_dim])/chunk_size[div_dim]+0.5);
    left->max[div_dim] = t->root->min[div_dim] + total_chunks*chunk_size[div_dim];
    right->min[div_dim] = left->max[div_dim];

    factors[div_dim] = left_chunks;
    last_block = _fast3tree_parallel_mapping(t, left, factors, factors_prod, chunk_size, mapping, div_dim, last_block);
    right = t->root + left_index + 1;
    factors[div_dim] = right_chunks;
    last_block = _fast3tree_parallel_mapping(t, right, factors, factors_prod, chunk_size, mapping, div_dim, last_block);
    factors[div_dim] = prev_factor;
    return last_block;
  }

  n = 0;
  //No divisions necessary; just need to figure out mapping.
  for (i=0; i<FAST3TREE_DIM; i++)
    n += factors_prod[i]*((int16_t)((0.5*(node->min[i]+node->max[i])-t->root->min[i])/chunk_size[i]));
  mapping[n] = last_block;
  return(last_block+1);
}

#undef _fast3tree_split_nodes_parallel
#define _fast3tree_split_nodes_parallel _F3TN(FAST3TREE_PREFIX,_fast3tree_split_nodes_parallel)
void _fast3tree_split_nodes_parallel(struct fast3tree *t, struct tree3_node *node) {
  struct tree3_node *null_ptr = NULL;
  if (node->num_points <= POINTS_PER_LEAF) {
    node->div_dim = -1;
    return;
  }
  int64_t cpus = 0;
#pragma omp parallel
  cpus = omp_get_num_threads();
  
  int64_t real_cpus = cpus;
  while (cpus*FAST3TREE_PARALLEL_CHUNK_SIZE < node->num_points) cpus *= 2;

  int64_t factors[FAST3TREE_DIM], factors_prod[FAST3TREE_DIM];
  _fast3tree_check_init(mapping, int64_t, cpus);
  float chunk_size[FAST3TREE_DIM];
  int64_t i, n=0;

#ifdef FAST3TREE_PRINT_TIMING
  int64_t t0 = mach_absolute_time();
#endif /* FAST3TREE_PRINT_TIMING */

  if (FAST3TREE_DIM==1) { factors[0] = cpus; }
  else {
    //Split along dimensions according to factorization of # of CPUs
    for (i = ceil(fabs(pow(cpus,1.0/(double)FAST3TREE_DIM))); i>0&&n<(FAST3TREE_DIM-1); i--)
      if ((cpus % i)==0) {
	cpus /= i;
	factors[n++] = i;
	i = ceil(pow(fabs(cpus),1.0/((double)FAST3TREE_DIM-(double)n)))+1;
      }
    factors[FAST3TREE_DIM-1] = cpus;
  }

  node = t->root;
  for (i=0; i<FAST3TREE_DIM; i++)
    chunk_size[i] = 1.0001*((node->max[i]-node->min[i])/(double)factors[i]);
  cpus = 1;
  for (i=0; i<FAST3TREE_DIM; i++) { factors_prod[i] = cpus; cpus*=factors[i]; }
  assert(cpus > 0);

  assert(cpus == _fast3tree_parallel_mapping(t, node, factors, factors_prod, chunk_size, mapping, -1, 0));
  _fast3tree_timing("Prep time");
  _fast3tree_check_init(seg, int16_t, node->num_points);
  _fast3tree_check_init(counts, int64_t, cpus*real_cpus);
  memset(counts, 0, sizeof(int64_t)*cpus*real_cpus);
  _fast3tree_check_init(indices, int64_t, cpus);

  //Store exact floating-point representation of minimums
  float *mins[FAST3TREE_DIM];
  for (i=0; i<FAST3TREE_DIM; i++) {
    mins[i] = _fast3tree_check_realloc(NULL, sizeof(float)*(factors[i]+1), "mins");
    for (n=0; n<factors[i]+1; n++) mins[i][n] = node->min[i] + n*chunk_size[i];
  }

  //Calculate volume partitions for each point and pre-sort.
#pragma omp parallel default(none) shared(node,seg,counts,chunk_size,factors_prod,mins,mapping,cpus,real_cpus)
  {
    int64_t seg_i, n, i, to_process = 0;
    int64_t tid = omp_get_thread_num();
    int64_t block_size = ceil((double)node->num_points/(double)(real_cpus*cpus));
    _fast3tree_check_init(my_ind, int64_t, cpus);
    _fast3tree_check_init(my_max_ind, int64_t, cpus);
    _fast3tree_check_init(my_counts, int64_t, cpus);
    memset(my_counts, 0, sizeof(int64_t)*cpus);
    for (seg_i=0; seg_i<cpus; seg_i++) {
      my_ind[seg_i] = (seg_i*real_cpus+tid)*block_size;
      my_max_ind[seg_i] = my_ind[seg_i]+block_size;
      if (my_max_ind[seg_i]>node->num_points)
	my_max_ind[seg_i] = node->num_points;
      if (my_max_ind[seg_i]>my_ind[seg_i])
	to_process += my_max_ind[seg_i]-my_ind[seg_i];
    }
    for (seg_i=0; seg_i<cpus; seg_i++) {
      int64_t n_start = my_ind[seg_i];
      int64_t max_n = my_max_ind[seg_i];
      for (n=n_start; n<max_n; n++) {
	seg[n] = 0;
	for (i=0; i<FAST3TREE_DIM; i++) {
	  int64_t chunk = (node->points[n].pos[i]-node->min[i])/chunk_size[i];
	  if (node->points[n].pos[i] < mins[i][chunk]) chunk--;
	  else if (node->points[n].pos[i] > mins[i][chunk+1]) chunk++;
	  seg[n] += factors_prod[i]*chunk;
	}
	seg[n] = mapping[seg[n]];
	my_counts[seg[n]]++;
	if (seg[n]==seg_i) my_ind[seg_i]++;
	else if (my_ind[seg[n]] < my_max_ind[seg[n]]) {
	  //Swap:
	  FAST3TREE_TYPE tmp = node->points[my_ind[seg[n]]];
	  node->points[my_ind[seg[n]]] = node->points[n];
	  node->points[n] = tmp;
	  seg[my_ind[seg[n]]] = seg[n];
	  my_ind[seg[n]]++;
	  n--;
	}
	else my_ind[seg_i]++; //Can't do anything about it
      }
    }
    for (i=0; i<cpus; i++) to_process -= my_counts[i];
    assert(to_process == 0);

    free(my_ind);
    free(my_max_ind);
    memcpy(counts+cpus*tid, my_counts, sizeof(int64_t)*cpus);
    free(my_counts);
  }
  free(mapping);
  for (i=0; i<FAST3TREE_DIM; i++) free(mins[i]);

  int64_t total_points = 0;
  for (i=0; i<cpus; i++) {
    for (n=1; n<real_cpus; n++)
      counts[i] += counts[i+n*cpus];
    total_points += counts[i];
  }
  assert(total_points == node->num_points);
  _fast3tree_timing("Counting/segmenting/pre-sorting");

  //Resort remaining points; cheap, so single-thread for now
  indices[0]=0;
  for (i=1; i<cpus; i++) indices[i] = indices[i-1]+counts[i-1];
  for (n=0; n<node->num_points; n++) {
    int16_t ts = seg[n];
    if (ts < 0) continue;
    if (indices[ts] != n) {
      FAST3TREE_TYPE tmp = node->points[indices[ts]];
      node->points[indices[ts]] = node->points[n];
      node->points[n] = tmp;
      seg[n] = seg[indices[ts]];
      seg[indices[ts]] = -1;
      indices[ts]++;
      n--;
    }
    while (seg[indices[ts]]==ts) { seg[indices[ts]] = -1; indices[ts]++; }
  }
  free(seg);
  _fast3tree_timing("Resort");

  //Build trees in parallel
  t->root = _fast3tree_check_realloc(t->root, sizeof(struct tree3_node)*t->num_nodes, "Tree nodes");
  memset(indices, 0, sizeof(int64_t)*cpus);
  for (i=1; i<cpus; i++) indices[i] = indices[i-1]+counts[i-1];
  _fast3tree_check_init(tt, struct fast3tree *, cpus);
#pragma omp parallel for default(none) shared(tt,counts,indices,t,cpus) private(i)
  for (i=0; i<cpus; i++) {
    tt[i] = _fast3tree_check_realloc(NULL,sizeof(struct fast3tree), "Allocating fast3tree.");
    memset(tt[i], 0, sizeof(struct fast3tree));
    tt[i]->points = t->root->points + indices[i];
    tt[i]->num_points = counts[i];
    _fast3tree_build(tt[i], 1);
  }
  free(indices);

  _fast3tree_timing("Subtree build");

  //Copy new nodes to original tree
  int64_t node_index = 0;
  int64_t max_nodes = t->num_nodes;
  int64_t p_so_far = 0;
  int64_t total_nodes = max_nodes;
  for (i=0; i<cpus; i++) total_nodes += tt[i]->num_nodes-1;
  t->root = _fast3tree_check_realloc(t->root, sizeof(struct tree3_node)*(total_nodes), "Tree nodes");  //Memory won't be allocated until pages are actually written
  t->allocated_nodes = total_nodes;
  _fast3tree_timing("Tree alloc");

  _fast3tree_check_init(num_nodes, int64_t, (cpus+1));
  _fast3tree_check_init(div_dim, int16_t, cpus);
  num_nodes[0] = t->num_nodes;
  for (i=1; i<=cpus; i++) {
    num_nodes[i] = num_nodes[i-1]+tt[i-1]->num_nodes-1;
    div_dim[i-1] = tt[i-1]->root->div_dim;
  }

#pragma omp parallel for default(none) shared(cpus,num_nodes,t,tt,null_ptr) private(i,n)
  for (i=0; i<cpus; i++) {
    struct tree3_node *dest = t->root+num_nodes[i];
    struct tree3_node *src = tt[i]->root+1;
    int64_t max_n = tt[i]->num_nodes-1;
    for (n=0; n<max_n; n++) {
      dest[n] = src[n];
#define REBUILD(x) x = t->root + (x-null_ptr) + (num_nodes[i] - 1);
      REBUILD(dest[n].left);
      REBUILD(dest[n].right);
      REBUILD(dest[n].parent);
#undef REBUILD
    }
    fast3tree_free(tt+i);
  }
  free(tt);

  for (i=0; i<cpus; i++) {
    while (node_index < max_nodes && t->root[node_index].div_dim > -1)
      node_index++;
    assert(node_index < max_nodes);
    struct tree3_node *nd = t->root + node_index;
    nd->points = t->root->points + p_so_far;
    nd->num_points = counts[i];
    p_so_far += counts[i];
    if (num_nodes[i] < num_nodes[i+1]) {
      nd->left = null_ptr + num_nodes[i];
      nd->right = null_ptr + num_nodes[i] + 1;
      nd->div_dim = div_dim[i];
      t->root[num_nodes[i]].parent = 
	t->root[num_nodes[i]+1].parent = t->root + node_index;
      assert(nd->num_points == t->root[num_nodes[i]].num_points + t->root[num_nodes[i]+1].num_points);
    }
    node_index++;
  }
  t->num_nodes = num_nodes[cpus];
  free(div_dim);
  free(counts);
  free(num_nodes);

  //Fix node point counts and pointers
  for (i=max_nodes-1; i>=0; i--) {
#define REBUILD(x) x = t->root + (x - null_ptr);
    REBUILD(t->root[i].parent);
    REBUILD(t->root[i].left);
    REBUILD(t->root[i].right);
#undef REBUILD
    if (i>0 && t->root[i].num_points == t->root->num_points) {
      t->root[i].points = t->root[i].left->points;
      t->root[i].num_points = t->root[i].left->num_points + t->root[i].right->num_points;
    }
  }
  assert(t->root->div_dim < 0 || (t->root->num_points == t->root->left->num_points+t->root->right->num_points));

  _fast3tree_timing("Subtree merge");
}

#endif /*OPENMP*/


#undef _fast3tree_rebuild_pointers
#define _fast3tree_rebuild_pointers _F3TN(FAST3TREE_PREFIX,_fast3tree_rebuild_pointers)
void _fast3tree_rebuild_pointers(struct fast3tree *t) {
  int64_t i;
  struct tree3_node *nullptr = NULL;
  for (i=0; i<t->num_nodes; i++) {
#define REBUILD(x) x = t->root + (x - nullptr)
    REBUILD(t->root[i].left);
    REBUILD(t->root[i].right);
    REBUILD(t->root[i].parent);
#undef REBUILD
  }
}

#undef _fast3tree_build
#define _fast3tree_build _F3TN(FAST3TREE_PREFIX,_fast3tree_build)
void _fast3tree_build(struct fast3tree *t, int64_t within_parallel) {
  int64_t i, j;
  struct tree3_node *root;
  FAST3TREE_TYPE tmp;
  t->allocated_nodes = (3+t->num_points/(POINTS_PER_LEAF/2));
  t->root = _fast3tree_check_realloc(t->root, sizeof(struct tree3_node)*(t->allocated_nodes), "Tree nodes"); //Estimate memory load
  t->num_nodes = 1;
  int64_t use_parallel_preprocess = 0;
  int64_t skip_rebuild = 0;
#ifdef _OPENMP
  if (!within_parallel) use_parallel_preprocess = 1;
#endif /* _OPENMP */

  if (!use_parallel_preprocess) {
    //Get rid of NaNs / infs
    for (i=0; i<t->num_points; i++) {
      for (j=0; j<FAST3TREE_DIM; j++) if (!isfinite(t->points[i].pos[j])) break;
      if (j<FAST3TREE_DIM) {
	tmp = t->points[i];
	t->num_points--;
	t->points[i] = t->points[t->num_points];
	t->points[t->num_points] = tmp;
	i--;
      }
    }
  }

  root = t->root;
  memset(root, 0, sizeof(struct tree3_node));
  root->num_points = t->num_points;
  root->points = t->points;
  root->div_dim = -1;
  if (t->num_points) {
    if (!use_parallel_preprocess) _fast3tree_find_minmax(root);
#ifdef _OPENMP
    else _fast3tree_preprocess_parallel(t, root);
#endif /* _OPENMP */
  }
  for (j=0; j<FAST3TREE_DIM; j++) assert(isfinite(root->min[j]));
  for (j=0; j<FAST3TREE_DIM; j++) assert(isfinite(root->max[j]));

#ifdef _OPENMP
  if (root->num_points > FAST3TREE_MIN_PARALLEL && !within_parallel) {
    _fast3tree_split_nodes_parallel(t, root);
    skip_rebuild = 1;
  }
  else
#endif /* _OPENMP */
  if (root->num_points > POINTS_PER_LEAF)
    _fast3tree_split_node(t, root);

  if (!skip_rebuild) {
    t->root = _fast3tree_check_realloc(t->root, sizeof(struct tree3_node)*(t->num_nodes), "Tree nodes");
    t->allocated_nodes = t->num_nodes;
    if (!within_parallel)
      _fast3tree_rebuild_pointers(t);
  }
}

#undef _fast3tree_maxmin_rebuild
#define _fast3tree_maxmin_rebuild _F3TN(FAST3TREE_PREFIX,_fast3tree_maxmin_rebuild)
void _fast3tree_maxmin_rebuild(struct tree3_node *n) {
  int i;
  if (n->div_dim < 0 && n->num_points) {
    _fast3tree_find_minmax(n);
    return;
  }
  _fast3tree_maxmin_rebuild(n->left);
  _fast3tree_maxmin_rebuild(n->right);
  memcpy(n->min, n->left->min, sizeof(float)*FAST3TREE_DIM);
  memcpy(n->max, n->right->max, sizeof(float)*FAST3TREE_DIM);
  for (i=0; i<FAST3TREE_DIM; i++) {
    if (n->min[i] > n->right->min[i]) n->min[i] = n->right->min[i];
    if (n->max[i] < n->left->max[i]) n->max[i] = n->left->max[i];
  }
}


#undef _fast3tree_check_realloc
#define _fast3tree_check_realloc _F3TN(FAST3TREE_PREFIX,_fast3tree_check_realloc)
void *_fast3tree_check_realloc(void *ptr, size_t size, char *reason) {
  void *res = realloc(ptr, size);
  if ((res == NULL) && (size > 0)) {
    fprintf(stderr, "[Error] Failed to allocate memory (%s)!\n", reason);
    exit(1);
  }
  return res;
}


#undef _fast3tree_set_minmax
#define _fast3tree_set_minmax _F3TN(FAST3TREE_PREFIX,fast3tree_set_minmax)
void _fast3tree_set_minmax(struct fast3tree *t, float min, float max) {
  int i;
  for (i=0; i<FAST3TREE_DIM; i++) {
    t->root->min[i] = min;
    t->root->max[i] = max;
  }
}


#undef _fast3tree_find_next_closest_dist
#define _fast3tree_find_next_closest_dist _F3TN(FAST3TREE_PREFIX,_fast3tree_find_next_closest_dist)
float _fast3tree_find_next_closest_dist(const struct tree3_node *n, const float c[FAST3TREE_DIM], float r, const struct tree3_node *o_n, int64_t *counts) {
  int64_t i,j;
  float r2, dist, dx, *pos;

  if (_fast3tree_box_not_intersect_sphere(n,c,r) || (o_n==n)) return r;
  if (n->div_dim < 0) { /* Leaf node */
    r2 = r*r;
    for (i=0; i<n->num_points; i++) {
      j = dist = 0;
      pos = n->points[i].pos;
      for (; j<FAST3TREE_DIM; j++) {
	dx = c[j]-pos[j];
	dist += dx*dx;
      }
      if (dist < r2) r2 = dist;
    }
    *counts = (*counts)+1;
    return sqrt(r2);
  }
  struct tree3_node *n1=n->left, *n2=n->right;
  if (c[n->div_dim] > 0.5*(n->min[n->div_dim]+n->max[n->div_dim])) {
    n1 = n->right;
    n2 = n->left;
  }
  r = _fast3tree_find_next_closest_dist(n1, c, r, o_n, counts);
  return _fast3tree_find_next_closest_dist(n2, c, r, o_n, counts);
}


#undef fast3tree_find_next_closest_distance
#define fast3tree_find_next_closest_distance _F3TN(FAST3TREE_PREFIX,fast3tree_find_next_closest_distance)
float fast3tree_find_next_closest_distance(struct fast3tree *t, struct fast3tree_results *res, float c[FAST3TREE_DIM]) {
  int64_t i=0, j;
  float dist = 0, dx, min_dist = 0;
  struct tree3_node *nd = t->root;
  
  while (nd->div_dim >= 0) {
    if (c[nd->div_dim] <= (nd->left->max[nd->div_dim])) { nd = nd->left; }
    else { nd = nd->right; }
  }

  while (nd != t->root && (nd->min[nd->parent->div_dim] == nd->max[nd->parent->div_dim])) nd = nd->parent;

  min_dist = 0;
  for (j=0; j<FAST3TREE_DIM; j++) {
    dx = nd->max[j]-nd->min[j];
    min_dist += dx*dx;
  }

  for (i=0; i<nd->num_points; i++) {
    dist = 0;
    for (j=0; j<FAST3TREE_DIM; j++) {
      dx = c[j] - nd->points[i].pos[j];
      dist += dx*dx;
    }
    if (dist && (dist < min_dist)) min_dist = dist;
  }
  int64_t counts = 0;
  min_dist = _fast3tree_find_next_closest_dist(t->root,c,sqrt(min_dist), nd, &counts);
  return min_dist;
}

#undef float

#endif /* _FAST3TREE_C_ */
