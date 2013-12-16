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

#ifndef CHARGEASSIGN_H
#define CHARGEASSIGN_H

#include "types.h"
#include "p3m-common.h"

#include "interpol.h"



void assign_charge(system_t *, parameters_t *, data_t *, int);
void assign_forces(FLOAT_TYPE, system_t *, parameters_t *, data_t *, forces_t *, int);
void assign_forces_real(FLOAT_TYPE force_prefac, system_t *s, parameters_t *p, data_t *d, forces_t *f);
void assign_forces_ad_real(double force_prefac, system_t *s, parameters_t *p, data_t *d, forces_t *f);

void assign_charge_real(system_t *s, parameters_t *p, data_t *d);
void assign_forces_real(FLOAT_TYPE force_prefac, system_t *s, parameters_t *p, data_t *d, forces_t *f);

void assign_forces_interlacing(FLOAT_TYPE force_prefac, system_t *s, parameters_t *p, data_t *d, forces_t *f);

void assign_forces_ad(double force_prefac, system_t* s, parameters_t* p, data_t* d, forces_t *, int ii);
void assign_charge_and_derivatives(system_t *, parameters_t *, data_t *, int);
void assign_charge_and_derivatives_real(system_t *s, parameters_t *p, data_t *d);

void assign_charge_q2(system_t *s, parameters_t *p, FLOAT_TYPE *Qmesh, int mesh, interpolation_t *inter);

void assign_charge_nocf(system_t *s, parameters_t *p, FLOAT_TYPE *Qmesh, int mesh, interpolation_t *inter);

void collect_rms_nocf(system_t *s, parameters_t *p, FLOAT_TYPE *Qmesh, FLOAT_TYPE *rms, int mesh, interpolation_t *inter);

#ifdef CA_DEBUG
#define CA_TRACE(A) A
#else
#define CA_TRACE
#endif

#endif
