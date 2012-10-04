#include <assert.h>
#include <math.h>

#include "domain-decomposition.h"
#include "common.h"

static void init_neighbors(domain_decomposition_t *d, cell_t *c) {
  celllist_t *i, *last=NULL, *buf;
  int n[3], m[3], ind=0;
  
  buf = Init_array( 26, sizeof(celllist_t));
  c->neighbors = buf;

  for(n[0]=(c->coords[0]-1);n[0]<=(c->coords[0]+1);n[0]++)
    for(n[1]=(c->coords[1]-1);n[1]<=(c->coords[1]+1);n[1]++)
      for(n[2]=(c->coords[2]-1);n[2]<=(c->coords[2]+1);n[2]++) {
	if((n[0]==c->coords[0]) &&
	   (n[1]==c->coords[1]) &&
	   (n[2]==c->coords[2]))
	  continue;

	for(int i = 0;i<3;i++) {
	  m[i] = n[i];
	  if(m[i] < 0)
	    m[i] += d->cells_per_direction;
	  if(m[i] >= d->cells_per_direction )
	    m[i] -= d->cells_per_direction;
	}
	if((m[0]==c->coords[0]) &&
	   (m[1]==c->coords[1]) &&
	   (m[2]==c->coords[2]))
	  continue;

	i = buf + ind;
	i->prev = (struct celllist_t *) last;
	if(last != NULL)
	  last->next = (struct celllist_t *)i;
	i->c = (struct cell_t *)  &(d->cells[d->cells_per_direction * d->cells_per_direction * m[0] +
			  d->cells_per_direction * m[1] + m[2]]);
	last = i;
	ind++;
      }

  last->next = NULL;
}

void Free_cell(cell_t *c) {
  assert( c != NULL);
  celllist_t *n = c->neighbors;
  while(n != NULL) {
    n = n->next;
    if(n->prev != NULL)
      FFTW_FREE(n->prev);
  }
  Resize_buffered_list(c->__q, 0);
  FFTW_FREE(c->__q);

  Resize_buffered_list(c->__ids, 0);
  FFTW_FREE(c->__ids);

  Resize_bvector_array( c->p, 0);
  FFTW_FREE(c->p);

  FFTW_FREE(c);
}

void Free_dd( domain_decomposition_t *dd) {
  assert( dd != NULL );

  for(int i=0;i<dd->total_cells;i++) {
    Free_cell(&dd->cells[i]);
  }

  FFTW_FREE(dd->cells);
  FFTW_FREE(dd);
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
	c->p = Init_bvector_array(0);
	c->__q = Init_buffered_list(0);
	c->q = c->__q->data;
	c->__ids = Init_buffered_list(0);
	c->ids = c->__ids->data;
	c->neighbors = NULL;
	c->n_particles = 0;
	for(int j=0;j<3;j++)
	  c->coords[j] = n[j];

	init_neighbors(d, c);
      }
  return d;
}

static inline void add_particle_to_cell( cell_t *c, int id, FLOAT_TYPE pos[3], FLOAT_TYPE q) {
  assert( c != 0 );
  int old_size = c->n_particles;
  int new_size = old_size+1;

  Resize_bvector_array( c->p, new_size );
  Resize_buffered_list(c->__q, new_size*sizeof(FLOAT_TYPE));
  c->q = c->__q->data;
  Resize_buffered_list(c->__ids, new_size*sizeof(int));
  c->ids = c->__ids->data;

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
