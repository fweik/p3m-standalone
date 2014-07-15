/**    Copyright (C) 2011,2012,2013 Florian Weik <fweik@icp.uni-stuttgart.de>

       This program is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

       This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this program.  If not, see <http://www.gnu.org/licenses/>. **/

#ifndef DOMAIN_DECOMPOSITION_H
#define DOMAIN_DECOMPOSITION_H

#include "types.h"

#define ILIST_STEP 100

typedef struct celllist_t celllist_t;

typedef struct {
  int real_size;
  int used_size;
  int *data;
} intlist_t;

typedef struct {
  int n_particles;
  bvector_array_t *p;
  bvector_array_t *v;
  FLOAT_TYPE *q;
  buffered_list_t *__q;
  int *ids;
  buffered_list_t *__ids;
  celllist_t *neighbors;
  int coords[3];
} cell_t;

struct celllist_t {
  cell_t *c;
  celllist_t *next;
  celllist_t *prev;
};

typedef struct {
  int cells_per_direction;
  int total_cells;
  cell_t *cells;
  FLOAT_TYPE h;
} domain_decomposition_t;

domain_decomposition_t *Init_dd( int cells_per_direction, FLOAT_TYPE box );
void Free_dd( domain_decomposition_t *dd);
cell_t *add_particle( domain_decomposition_t *d, int id, FLOAT_TYPE pos[3], FLOAT_TYPE q);
void add_system( domain_decomposition_t *d, system_t *s);

#endif
