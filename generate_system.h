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

#ifndef GENERATE_SYSTEM_H
#define GENERATE_SYSTEM_H

#include "common.h"

enum { SYSTEM_RANDOM = 0,
       SYSTEM_INNER_BOX = 1,
       SYSTEM_MADELUNG = 2,
       SYSTEM_SEPARATED_DIPOLE = 3,
       SYSTEM_GAUSSIAN = 4,
       SYSTEM_SLAB = 5
     };

system_t *generate_system( int, int, FLOAT_TYPE, FLOAT_TYPE);

#endif
