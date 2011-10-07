#ifndef COMMON_H
#define COMMON_H

#include <stdlib.h>

#include "types.h"

static inline FLOAT_TYPE my_power( FLOAT_TYPE base, int exponent ) {
  FLOAT_TYPE ret = 1.0;

  while( exponent-- > 0 ) {
    ret *= base;
  }

  return ret;
}



system_t *Init_system(int);
void *Init_array(int, size_t);vector_array_t *Init_vector_array(int);
forces_t *Init_forces(int);
void Free_system(system_t *);
void Free_forces(forces_t *);
void Free_vector_array(vector_array_t *);

void Calculate_forces ( const method_t *, system_t *, parameters_t *, data_t *, forces_t * );
void Calculate_reference_forces ( system_t *, parameters_t * );

#endif
