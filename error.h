#ifndef ERROR_H
#define ERROR_H

#include "common.h"

typedef struct {
  FLOAT_TYPE f;
  FLOAT_TYPE f_k;
  FLOAT_TYPE f_r;
  FLOAT_TYPE f_v[3];
} error_t;

error_t Calculate_errors(system_t *, forces_t *);

#endif