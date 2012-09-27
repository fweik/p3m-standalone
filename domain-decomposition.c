#include <assert.h>
#include <math.h>

#include "domain-decomposition.h"
#include "common.h"

void init_neighbors(domain_decomposition_t *d, cell_t *c) {
  celllist_t *i, *last=NULL;
  int n[3];

  for(n[0]=(c->coords[0]-1);n[0]<=(c->coords[0]+1);n[0]++)
    for(n[1]=(c->coords[1]-1);n[1]<=(c->coords[1]+1);n[1]++)
      for(n[2]=(c->coords[2]-1);n[2]<=(c->coords[2]+1);n[2]++) {
	if((n[0]==c->coords[0]) &&
	   (n[1]==c->coords[1]) &&
	   (n[2]==c->coords[2]))
	  continue;

	i = Init_array( 1, sizeof(celllist_t));
	i->prev = (struct celllist_t *) last;
	last->next = (struct celllist_t *)i;
	i->c = &(d->cells[d->cells_per_direction * d->cells_per_direction * n[0] +
			  d->cells_per_direction * n[1] + n[2]]);
	last = i;
      }

  last->next = NULL;
}

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
	cell_t *c = &(d->cells[ind++]);
	c->p = NULL;
	c->q = NULL;
	c->ids = NULL;
	c->neighbors = NULL;
	c->n_particles = 0;
	for(int j=0;j<3;j++)
	  c->coords[j] = n[j];

	init_neighbors(d, c);
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
    c->p->fields[i][old_size] = pos[i];

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
