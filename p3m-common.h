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

#ifndef P3M_COMMON_H
#define P3M_COMMON_H

#include "types.h"
#include "wtime.h"

extern int P3M_BRILLOUIN_TUNING;
extern int P3M_BRILLOUIN;

#define r_ind(A,B,C) ((A)*d->mesh*d->mesh + (B)*d->mesh + (C))
#define c_ind(A,B,C) (2*d->mesh*d->mesh*(A)+2*d->mesh*(B)+2*(C))

FLOAT_TYPE sinc(FLOAT_TYPE);
FLOAT_TYPE analytic_cotangent_sum(int n, FLOAT_TYPE mesh_i, int cao);

void Init_differential_operator( data_t * );
void Init_nshift(data_t *);
data_t *Init_data(const method_t *, system_t *s, parameters_t *); 
void Free_data(data_t *);

FLOAT_TYPE C_ewald(int nx, int ny, int nz, system_t *s, parameters_t *p);
FLOAT_TYPE C_ewald_dip(int nx, int ny, int nz, system_t *s, parameters_t *p);
FLOAT_TYPE C_ewald_water(int nx, int ny, int nz, system_t *s, parameters_t *p);
FLOAT_TYPE Generic_error_estimate(R3_to_R A, R3_to_R B, R3_to_R C, system_t *s, parameters_t *p, data_t *d);

FLOAT_TYPE Generic_error_estimate_inhomo(system_t *s, parameters_t *p, int uniform, int mesh, int cao, int mc, char *out_file);

FLOAT_TYPE A_const(int nx, int ny, int nz, system_t *s, parameters_t *p);
FLOAT_TYPE B_const(int nx, int ny, int nz, system_t *s, parameters_t *p);

FLOAT_TYPE *Error_map(system_t *s, forces_t *f, forces_t *f_ref, int mesh, int cao);

#define TIMING_START_C TIMING_START(t_c)
#define TIMING_START_F TIMING_START(t_f)
#define TIMING_START_G TIMING_START(t_g)

#define TIMING_STOP_C TIMING_STOP(t_c)
#define TIMING_STOP_F TIMING_STOP(t_f)
#define TIMING_STOP_G TIMING_STOP(t_g)

#define TIMING_START(A) if(p->tuning) d->runtime.A = wtime();
#define TIMING_STOP(A) if(p->tuning) d->runtime.A = wtime() - d->runtime.A;

#endif
