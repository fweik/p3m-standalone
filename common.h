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

#define __detailed_timings

#ifdef __detailed_timings
extern double t_charge_assignment[4];
extern double t_force_assignment[4];
extern double t_real_part[4];
extern double t_k_part[4];
extern double t_convolution[4];
extern double t_fft[4];
extern double timer;
#endif

system_t *Init_system(int);
void Free_system(system_t *);

forces_t *Init_forces(int);
void Free_forces(forces_t *);

void *Init_array(int, size_t);
void *Resize_array(void *a, size_t new_size, size_t old_size);

vector_array_t *Init_vector_array(int);
void Resize_vector_array(vector_array_t *d, int new_size);
void Free_vector_array(vector_array_t *);

void Calculate_forces ( const method_t *, system_t *, parameters_t *, data_t *, forces_t * );
FLOAT_TYPE Calculate_reference_forces ( system_t *, parameters_t * );

FLOAT_TYPE Min_distance( system_t *s);

#endif
