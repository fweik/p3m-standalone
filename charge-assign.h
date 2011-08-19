#ifndef CHARGEASSIGN_H
#define CHARGEASSIGN_H

#include "p3m.h"

#include "interpol.h"

void assign_charge(system_t *, p3m_parameters_t *, p3m_data_t *, int);
void assign_forces(FLOAT_TYPE, system_t *, p3m_parameters_t *, p3m_data_t *, int);

void assign_forces_ad(system_t *s, p3m_parameters_t *p, p3m_data_t *d, int ii);
void assign_charge_and_derivatives(system_t *s, p3m_parameters_t *p, p3m_data_t *d, int ii);

#ifdef CA_DEBUG
#define CA_TRACE(A) A
#else
#define CA_TRACE
#endif

#endif
