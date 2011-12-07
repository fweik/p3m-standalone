#ifndef TUNING_H
#define TUNING_H

#include "types.h"

// @TODO: Make this more elegant

// Mesh sizes suitable for FFT, these are the 13-smooth numbers, excluding those which have prime factors 11 and 13.

#define MESH_MIN 7
#define MESH_MAX 256
#define MESH_STEP 4

#define ALPHA_STEP_MIN 1e-4
#define ALPHA_STEP 0.4

#define CAO_MIN 2
#define CAO_MAX 7

int Tune( const method_t *, system_t *, parameters_t *, FLOAT_TYPE );

#endif
