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

interpolation_t *Init_interpolation(int, int);
void Free_interpolation(interpolation_t *i);

#endif
