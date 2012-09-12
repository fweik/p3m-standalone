#ifndef DOMAIN_DECOMPOSITION_H
#define DOMAIN_DECOMPOSITION_H

#include "types.h"

typedef struct {
  int n_particles;
  vector_array_t *p;
  FLOAT_TYPE *q;
  int *ids;
  struct cell_t *neighbors;
  int coords[3];
} cell_t;

typedef struct {
  int cells_per_direction;
  int total_cells;
  cell_t *cells;
  FLOAT_TYPE h;
} domain_decomposition_t;

#endif
