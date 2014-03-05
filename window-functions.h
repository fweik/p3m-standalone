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
