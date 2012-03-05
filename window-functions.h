#ifndef WINDOW_FUNCTIONS_H
#define WINDOW_FUNCTIONS_H

#include "types.h"

/** Computes the  assignment function of for the \a i'th degree
    at value \a x. */
FLOAT_TYPE caf_bspline_k(int i, FLOAT_TYPE d);
FLOAT_TYPE caf_bspline(int i, FLOAT_TYPE x, int cao_value);
FLOAT_TYPE caf_bspline_d(int i, FLOAT_TYPE x, int cao_calue);

FLOAT_TYPE caf_kaiserbessel_k(int i, FLOAT_TYPE d);
FLOAT_TYPE caf_kaiserbessel(int i, FLOAT_TYPE x, int cao);

#endif
