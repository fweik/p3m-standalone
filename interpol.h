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
