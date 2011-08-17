#pragma once

#ifndef P3M_COMMON_H
#define P3M_COMMON_H

#include "p3m.h"

#define P3M_BRILLOUIN_TUNING 0
#define P3M_BRILLOUIN 0

#define PI 3.14159265358979323846264

#define r_ind(A,B,C) ((A)*Mesh*Mesh + (B)*Mesh + (C))
#define c_ind(A,B,C) (2*Mesh*Mesh*(A)+2*Mesh*(B)+2*(C))

FLOAT_TYPE sinc(FLOAT_TYPE);

void Differenzenoperator_berechnen(void);
void nshift_ausrechnen(void);
#endif
