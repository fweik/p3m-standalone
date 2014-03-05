/**    Copyright (C) 2011,2012,2013 Florian Weik <fweik@icp.uni-stuttgart.de>

       This program is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

       This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this program.  If not, see <http://www.gnu.org/licenses/>. **/

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
void Free_array(void *a);

buffered_list_t *Init_buffered_list(size_t size);
void Resize_buffered_list(buffered_list_t *l, size_t new_size);

vector_array_t *Init_vector_array(int);
void Resize_vector_array(vector_array_t *d, int new_size);
void Free_vector_array(vector_array_t *);

bvector_array_t *Init_bvector_array(int n);
void Resize_bvector_array(bvector_array_t *d, int new_size);

void Calculate_forces ( const method_t *, system_t *, parameters_t *, data_t *, forces_t * );
FLOAT_TYPE Calculate_reference_forces ( system_t *, parameters_t * );

FLOAT_TYPE Min_distance( system_t *s);

#endif
