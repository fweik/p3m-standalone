#ifndef TUNING_H
#define TUNING_H

#include "types.h"

// @TODO: Make this more elegant

#define CAO_MIN 2
#define CAO_MAX 7

#define N_TUNING_SAMPLES 100

typedef struct {
  double avg;
  double sgm;
  int n;
} timing_t;

timing_t Tune( const method_t *, system_t *, parameters_t *, FLOAT_TYPE );

#endif
