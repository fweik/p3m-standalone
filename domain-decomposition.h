#ifndef DOMAIN_DECOMPOSITION_H
#define DOMAIN_DECOMPOSITION_H

#include "types.h"

typedef struct {
  struct cell_t *c;
  struct celllist_t *next;
  struct celllist_t *prev;
  int wraped[3];
} celllist_t;

typedef struct {
  int n_particles;
  bvector_array_t *p;
  bvector_array_t *v;
  FLOAT_TYPE *q;
  buffered_list_t *__q;
  int *ids;
  buffered_list_t *__ids;
  struct celllist_t *neighbors;
  int coords[3];
} cell_t;

typedef struct {
  int cells_per_direction;
  int total_cells;
  cell_t *cells;
  FLOAT_TYPE h;
} domain_decomposition_t;

domain_decomposition_t *Init_dd( int cells_per_direction, FLOAT_TYPE box );
cell_t *add_particle( domain_decomposition_t *d, int id, FLOAT_TYPE pos[3], FLOAT_TYPE q);

#endif
