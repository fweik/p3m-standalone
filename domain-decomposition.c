#include <assert.h>
#include <math.h>

#include "domain-decomposition.h"
#include "common.h"

domain_decomposition_t *Init_dd( int cells_per_direction, FLOAT_TYPE box ) {
  domain_decomposition_t *d = Init_array( 1, sizeof(domain_decomposition_t) );
  int n[3], ind=0;

  assert(cells_per_direction != 0);

  d->h = box / cells_per_direction;
  
  d->total_cells = cells_per_direction * cells_per_direction * cells_per_direction;

  d->cells_per_direction = cells_per_direction;

  d->cells = Init_array( d->total_cells, sizeof(cell_t));

  for(n[0]=0;n[0]<cells_per_direction;n[0]++)
    for(n[1]=0;n[1]<cells_per_direction;n[1]++)
      for(n[2]=0;n[2]<cells_per_direction;n[2]++) {
	cell_t *c = &(d->cells[ind]);
	c->p = NULL;
	c->q = NULL;
	c->ids = NULL;
	c->n_particles = 0;
	for(int j=0;j<3;j++)
	  c->coords[j] = n[j];
      }
  return d;
}

void add_particle_to_cell( cell_t *c, int id, FLOAT_TYPE pos[3], FLOAT_TYPE q) {
  assert( c != 0 );
  int old_size = c->n_particles;
  int new_size = old_size+1;
   
  Resize_vector_array( c->p, new_size );
  c->q = Resize_array( c->q, new_size*sizeof(FLOAT_TYPE), old_size*sizeof(FLOAT_TYPE));
  c->ids = Resize_array( c->ids, new_size*sizeof(int), old_size*sizeof(int));

  for(int i = 0; i<3; i++)
    c->p->fields[i][old_size] = p[i];

  c->q[old_size] = q;
  c->ids[old_size] = id;

  c->n_particles++;     
}

cell_t *add_particle( domain_decomposition_t *d, int id, FLOAT_TYPE pos[3], FLOAT_TYPE q) {
  int n[3], ind;
  cell_t *c;

  assert(d != NULL);

  for(int i=0;i<3;i++)
    n[i] = (int)FLOOR(pos[i] / d->h);

  ind = d->cells_per_direction * d->cells_per_direction * n[0] +
    d->cells_per_direction * n[1] + n[2];

  if(ind >= d->total_cells)
    return NULL;

  c = &(d->cells[ind]);

  add_particle_to_cell( c, id, pos, q);

  return c;
}
