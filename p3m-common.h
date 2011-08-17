#ifndef P3M_COMMON_H
#define P3M_COMMON_H

#include "p3m.h"

#define P3M_BRILLOUIN_TUNING 0
#define P3M_BRILLOUIN 0

#define PI 3.14159265358979323846264

#define r_ind(A,B,C) ((A)*d->mesh*d->mesh + (B)*d->mesh + (C))
#define c_ind(A,B,C) (2*d->mesh*d->mesh*(A)+2*d->mesh*(B)+2*(C))

FLOAT_TYPE sinc(FLOAT_TYPE);

void Differenzenoperator_berechnen(p3m_parameters_t *, p3m_data_t *);
void nshift_ausrechnen(p3m_data_t *);
#endif
