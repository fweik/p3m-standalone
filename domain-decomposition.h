#ifndef DOMAIN_DECOMPOSITION_H
#define DOMAIN_DECOMPOSITION_H

#include "types.h"

typedef struct {
  int n_particles;
  vector_array_t *p;
  FLOAT_TYPE *q;
  int *ids;
  struct celllist_t *neighbors;
  int coords[3];
} cell_t;

typedef struct {
  cell_t *c;
  struct celllist_t *next;
  struct celllist_t *prev;
  int pbc[3];
} celllist_t;

typedef struct {
  int cells_per_direction;
  int total_cells;
  cell_t *cells;
  FLOAT_TYPE h;
} domain_decomposition_t;

#endif
