#include "generate_system.h"

#include "common.h"

#include "math.h"

#include "stdlib.h"

system_t *generate_random_system(int size, FLOAT_TYPE box, FLOAT_TYPE max_charge ) {
  int i,j;
  system_t *s;

  s = Init_system(size);
  s->length = box;
  s->q2 = 0.0;
  
  for(i=0;i<size;i++) {
    for(j=0;j<2;j++) {
      s->p->fields[j][i] = box*drand48();
    }
    s->q[i] = 1.0 - 2.0 * (i%2);
    s->q2 += s->q[i] * s->q[i];
  }
  return s;
}

system_t *generate_system( int form_factor, int size, FLOAT_TYPE box, FLOAT_TYPE max_charge ) {
  if( form_factor == FORM_FACTOR_RANDOM ) {
    return generate_random_system( size, box, max_charge );
  }
  return NULL;
}