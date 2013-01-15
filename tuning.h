#ifndef TUNING_H
#define TUNING_H

#include "types.h"

// @TODO: Make this more elegant

#define CAO_MIN 3
#define CAO_MAX 7

int Tune( const method_t *, system_t *, parameters_t *, FLOAT_TYPE );

#endif
