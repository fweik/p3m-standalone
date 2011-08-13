#ifndef CHARGEASSIGN_H
#define CHARGEASSIGN_H

#include "p3m.h"

void assign_charge(int, FLOAT_TYPE, FLOAT_TYPE*, FLOAT_TYPE*, int);
void assign_forces(FLOAT_TYPE, FLOAT_TYPE*, int, FLOAT_TYPE*, int);

void assign_forces_ad(double, FLOAT_TYPE **, int, FLOAT_TYPE *,int);
void assign_charge_and_derivatives(int, FLOAT_TYPE, FLOAT_TYPE *, FLOAT_TYPE *, int);
#endif
