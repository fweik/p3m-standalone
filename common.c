#include <stdlib.h>
#include <assert.h>

#include "common.h"

void Init_array(void *a, int size, size_t field_size) {
  assert(size > 0);
  assert(field_size > 0);

  a = malloc(size * field_size);

  assert(a != NULL);
}

void Init_vector_array(vector_array_t *v, int n) {
  int i;

  assert(v != NULL);
  assert(n > 0 );

  for(i=0;i<3;i++) {
    Init_array(v->fields[i], n, sizeof(FLOAT_TYPE));
  }

  v->x = v->fields[0];
  v->y = v->fields[1];
  v->z = v->fields[2];
}

void Init_system(system_t *s) {
  assert(s!= NULL);
  assert( s->nparticles > 0 );

  Init_vector_array(&(s->p), s->nparticles);
  Init_vector_array(&(s->f), s->nparticles);
  Init_vector_array(&(s->f_k), s->nparticles);
  Init_vector_array(&(s->f_r), s->nparticles);

  Init_vector_array(&(s->reference.f), s->nparticles);
  Init_vector_array(&(s->reference.f_k), s->nparticles);

  Init_array(s->q, s->nparticles, sizeof(FLOAT_TYPE));

  assert((s->q) != NULL);
}
