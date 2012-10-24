#ifndef P3M_AD_SELF_FORCES_H
#define P3M_AD_SELF_FORCES_H

#include "types.h"

void Init_self_forces( system_t *s, parameters_t *p, data_t *d );
void Substract_self_forces( system_t *s, parameters_t *p, data_t *d, forces_t *f );

#endif
