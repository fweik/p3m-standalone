#ifndef STATISTICS_H
#define STATISTICS_H
#include <math.h>
#include <string.h>

#include "types.h"

FLOAT_TYPE *radial_charge_distribution(FLOAT_TYPE r_min, FLOAT_TYPE r_max, int bins, system_t *s);
FLOAT_TYPE *rdf_fft( int N, FLOAT_TYPE *rdf);
#endif
