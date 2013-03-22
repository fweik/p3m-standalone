#ifndef STATISTICS_H
#define STATISTICS_H
#include <math.h>
#include <string.h>

#include "types.h"

FLOAT_TYPE *radial_charge_distribution(FLOAT_TYPE r_min, FLOAT_TYPE r_max, int bins, system_t *s);
FLOAT_TYPE *rdf_fft( int N, FLOAT_TYPE *rdf);
FLOAT_TYPE *low_pass_forward( int N, FLOAT_TYPE *data, FLOAT_TYPE alpha);
FLOAT_TYPE *low_pass_backward( int N, FLOAT_TYPE *data, FLOAT_TYPE alpha);
void rshif_array( int N, FLOAT_TYPE *data, int shift );
void radial_distribution_species(FLOAT_TYPE r_min, FLOAT_TYPE r_max, int bins, system_t *s);
FLOAT_TYPE *radial_distribution(FLOAT_TYPE r_min, FLOAT_TYPE r_max, int bins, system_t *s);
#endif
