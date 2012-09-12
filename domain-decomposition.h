#ifndef DOMAIN_DECOMPOSITION_H
#define DOMAIN_DECOMPOSITION_H

#include "types.h"

#define ILIST_STEP 100

typedef struct {
  int real_size;
  int used_size;
  int *data;
} intlist_t;

typedef struct {
  int n_particles;
  vector_array_t *p;
  FLOAT_TYPE *q;
  intlist_t ids;
  struct cell_t *neighbors;
  int coords[3];
} cell_t;

typedef struct {
  int cells_per_direction;
  int total_cells;
  cell_t *cells;
  FLOAT_TYPE h;
} domain_decomposition_t;

domain_decomposition_t *Init_dd( int cells_per_direction, FLOAT_TYPE box );
void add_system( domain_decomposition_t *d, system_t *s);

#endif
