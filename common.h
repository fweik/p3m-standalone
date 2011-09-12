#ifndef COMMON_H
#define COMMON_H

#include <stdlib.h>

#include "types.h"

double my_power( double, int );

system_t *Init_system(int);
void *Init_array(int, size_t);
vector_array_t *Init_vector_array(int);
forces_t *Init_forces(int);
void Free_system(system_t *);
void Free_forces(forces_t *);
void Free_vector_array(vector_array_t *);

void Calculate_forces ( const method_t *, system_t *, parameters_t *, data_t *, forces_t * );
void Calculate_reference_forces ( system_t *, parameters_t * );

#endif
