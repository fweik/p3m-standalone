#ifndef TUNING_H
#define TUNING_H

#include "types.h"

// @TODO: Make this more elegant

#define MESH_MIN 4
#define MESH_MAX 256

#define RCUT_STEP 0.02
#define RCUT_STEP_MIN 0.001

#define CAO_MIN 1
#define CAO_MAX 7

parameters_t *Tune( const method_t *, system_t *, FLOAT_TYPE );

#endif
