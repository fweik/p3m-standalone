#include <stdlib.h>
#include <assert.h>

#include <stdio.h>

#include "common.h"

void *Init_array(int size, size_t field_size) {
  void *a;
  
  assert(size >= 0);
  assert(field_size > 0);

  a = malloc(size * field_size);

  assert(a != NULL);
  return a;
}

void Init_vector_array(vector_array_t *v, int n) {
  int i;

  assert(v != NULL);
  assert(n >= 0 );
  
  v->fields = Init_array( 3, sizeof(FLOAT_TYPE *));

  
  assert( v->fields != NULL );
  
  for(i=0;i<3;i++) {
    v->fields[i] = Init_array( n, sizeof(FLOAT_TYPE));
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

  s->q = Init_array(s->nparticles, sizeof(FLOAT_TYPE));

  assert((s->q) != NULL);
}

void Free_vector_array(vector_array_t *v) {
  int i;
      
  if( v != NULL ) {

    for(i=0;i<3;i++) {
      if(v->fields[i] != NULL) {
        free(v->fields[i]);
        v->fields[i] = NULL;
      }
    }

    v->x = v->y = v->z = NULL;
    v->fields = NULL;
  }
}