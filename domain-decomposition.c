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
	c->p = Init_vector_array(0);
	c->q = NULL;
	c->ids.data = NULL;
	c->ids.real_size = 0;
	c->ids.used_size = 0;
	c->n_particles = 0;
	for(int j=0;j<3;j++)
	  c->coords[j] = n[j];
	ind++;
      }
  return d;
}

void add_particle_to_cell( cell_t *c, int id, FLOAT_TYPE pos[3], FLOAT_TYPE q) {
  assert( c != 0 );
  int old_size = c->n_particles;
  int new_size = old_size+1;
   
  Resize_vector_array( c->p, new_size );
  c->q = Resize_array( c->q, new_size*sizeof(FLOAT_TYPE), old_size*sizeof(FLOAT_TYPE));

  if(new_size > c->ids.real_size) {
    c->ids.data = Resize_array( c->ids.data, (old_size+ILIST_STEP)*sizeof(int), old_size*sizeof(int));
    c->ids.real_size += ILIST_STEP;
  }


  for(int i = 0; i<3; i++)
    c->p->fields[i][old_size] = pos[i];

  c->q[old_size] = q;

  c->ids.data[old_size] = id;

  c->n_particles++;     
}

cell_t *add_particle( domain_decomposition_t *d, int id, FLOAT_TYPE pos[3], FLOAT_TYPE q) {
  int n[3], ind;
  cell_t *c;

  assert(d != NULL);

  for(int i=0;i<3;i++)
    n[i] = (int)FLOOR(pos[i] / d->h);

  /* printf("particle %d pos (%lf %lf %lf) cell (%d %d %d) q %lf\n", id, */
  /* 	 pos[0], pos[1], pos[2], n[0], n[1], n[2], q); */

  ind = d->cells_per_direction * d->cells_per_direction * n[0] +
    d->cells_per_direction * n[1] + n[2];

  if(ind >= d->total_cells)
    return NULL;

  c = &(d->cells[ind]);

  add_particle_to_cell( c, id, pos, q);

  return c;
}

void add_system( domain_decomposition_t *d, system_t *s) {
  FLOAT_TYPE pos[3];
  for(int id=0;id<s->nparticles;id++) {
    for(int j=0;j<3;j++)
      pos[j] = s->p->fields[j][id];
    add_particle( d, id, pos, s->q[id]);
  }
}
