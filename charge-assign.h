#ifndef CHARGEASSIGN_H
#define CHARGEASSIGN_H

#include "types.h"
#include "p3m-common.h"

#include "interpol.h"



void assign_charge(system_t *, parameters_t *, data_t *, int);
void assign_forces(FLOAT_TYPE, system_t *, parameters_t *, data_t *, forces_t *, int);
void assign_forces_real(FLOAT_TYPE force_prefac, system_t *s, parameters_t *p, data_t *d, forces_t *f);

void assign_forces_ad(double force_prefac, system_t* s, parameters_t* p, data_t* d, forces_t *, int ii);
void assign_charge_and_derivatives(system_t *, parameters_t *, data_t *, int);

void assign_charge_q2(system_t *s, parameters_t *p, FLOAT_TYPE *Qmesh, int mesh, interpolation_t *inter);

void assign_charge_nocf(system_t *s, parameters_t *p, FLOAT_TYPE *Qmesh, int mesh, interpolation_t *inter);

void collect_rms_nocf(system_t *s, parameters_t *p, FLOAT_TYPE *Qmesh, FLOAT_TYPE *rms, int mesh, interpolation_t *inter);

#ifdef CA_DEBUG
#define CA_TRACE(A) A
#else
#define CA_TRACE
#endif

#endif
