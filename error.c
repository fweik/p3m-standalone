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

#include <math.h>
#include <string.h>

#include "error.h"

error_t Calculate_errors(system_t *system, forces_t *f ) {
    int i,j;
    error_t e;
    memset(&e, 0, sizeof(error_t));

    for (i=0; i<system->nparticles; i++) {
        for (j=0;j<3;j++) {
            e.f   += SQR(  f->f->fields[j][i] - system->reference->f->fields[j][i] );
        }
    }

    e.f = SQRT(e.f);
    return e;
}
