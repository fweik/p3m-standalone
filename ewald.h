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

#ifndef EWALD_H
#define EWALD_H

#include "types.h"

#ifndef r_ind
#define r_ind(A,B,C) ((A)*d->mesh*d->mesh + (B)*d->mesh + (C))
#endif

data_t *Ewald_init(system_t *, parameters_t *);
void Ewald_k_space(system_t *, parameters_t *, data_t *, forces_t *);
void Ewald_compute_influence_function(system_t *, parameters_t *, data_t *);
FLOAT_TYPE Ewald_compute_optimal_alpha(system_t *, parameters_t *);
FLOAT_TYPE Ewald_estimate_error(system_t *, parameters_t *);
FLOAT_TYPE Ewald_error_k( system_t *, parameters_t *);

FLOAT_TYPE Ewald_energy(system_t *, parameters_t *);

extern const method_t method_ewald;

#endif
