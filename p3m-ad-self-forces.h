#ifndef P3M_AD_SELF_FORCES_H
#define P3M_AD_SELF_FORCES_H

#include "types.h"

//#define P3M_AD_SELF_FORCES_DEBUG

#ifdef P3M_AD_SELF_FORCES_DEBUG
#define SF_TRACE(A); A;
#else
#define SF_TRACE(A)
#endif


void Init_self_forces( system_t *s, parameters_t *p, data_t *d );
void Substract_self_forces( system_t *s, parameters_t *p, data_t *d, forces_t *f );

#endif
