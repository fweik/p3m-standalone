/**    Copyright (C) 2011,2012,2013 Florian Weik <fweik@icp.uni-stuttgart.de>

       This program is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

       This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this program.  If not, see <http://www.gnu.org/licenses/>. **/

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
