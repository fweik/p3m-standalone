#ifndef INTERPOL_H
#define INTERPOL_H

#include "types.h"
#include "window-functions.h"

#define MaxInterpol (2*100096)

/* Interpolation types */
enum {
  INTERPOL_BSPLINE,
  INTERPOL_KAISER
};

typedef struct {
  int  id;
  FLOAT_TYPE (*U)(int, FLOAT_TYPE, int);
  FLOAT_TYPE (*U_d)(int, FLOAT_TYPE, int);
  FLOAT_TYPE (*U_hat)(int, FLOAT_TYPE);
} interpolation_function_t;

interpolation_t *Init_interpolation(int, int);
void Free_interpolation(interpolation_t *i);

#endif
