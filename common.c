#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "p3m.h"
#include "common.h"

void *Init_array(int size, size_t field_size) {
    void *a;

    assert(size >= 0);
    assert(field_size > 0);

    a = malloc(size * field_size);

    assert(a != NULL);
    return a;
}

vector_array_t *Init_vector_array(int n) {
    int i;
    vector_array_t *v;

    v = Init_array( 1, sizeof(vector_array_t));

    assert(n >= 0 );

    v->fields = Init_array( 3, sizeof(FLOAT_TYPE *));


    assert( v->fields != NULL );

    for (i=0;i<3;i++) {
        v->fields[i] = Init_array( n, sizeof(FLOAT_TYPE));
    }

    v->x = v->fields[0];
    v->y = v->fields[1];
    v->z = v->fields[2];

    return v;
}

system_t *Init_system(int n) {
    system_t *s;

    assert( n > 0 );

    s = Init_array( 1, sizeof(system_t));

    s-> nparticles = n;

    s-> p = Init_vector_array(s->nparticles);

    s->reference = Init_forces(s->nparticles);

    s->q = Init_array(s->nparticles, sizeof(FLOAT_TYPE));

    return s;
}

forces_t *Init_forces(int n) {
    forces_t *f;

    f = Init_array( 1, sizeof(forces_t));

    f->f = Init_vector_array(n);
    f->f_k = Init_vector_array(n);
    f->f_r = Init_vector_array(n);
    
    return f;
}

void Free_system(system_t *s) {
  if( s == NULL)
    return;
  
  Free_vector_array(s->p);
  
  Free_forces(s->reference);
  
  free(s->q);
  
  free(s);
  
}

void Free_forces(forces_t *f) {
  if( f == NULL )
    return;
  
  Free_vector_array(f->f);
  Free_vector_array(f->f_k);
  Free_vector_array(f->f_r);
 
  free(f);
}

void Free_vector_array(vector_array_t *v) {
    int i;

    if ( v != NULL ) {
        for (i=0;i<3;i++) {
            if (v->fields[i] != NULL) {
                free(v->fields[i]);
                v->fields[i] = NULL;
            }
        }

        v->x = v->y = v->z = NULL;
        v->fields = NULL;
    }
    free(v);
}

