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

#ifndef TUNING_H
#define TUNING_H

#include "types.h"

// @TODO: Make this more elegant

#define CAO_MIN 2
#define CAO_MAX 7

#define N_TUNING_SAMPLES 1

typedef struct {
  double avg;
  double sgm;
  int n;
} timing_t;

timing_t Tune( const method_t *, system_t *, parameters_t *, FLOAT_TYPE );

#endif
