#ifndef P3M_IK_REAL_H
#define P3M_IK_REAL_H

#include "types.h"

void Influence_function_berechnen_ik_r(system_t*, parameters_t*, data_t*);
void P3M_ik_r(system_t *, parameters_t *, data_t *, forces_t *);
data_t *Init_ik_r(system_t*, parameters_t*);
FLOAT_TYPE Error_ik_r( system_t *, parameters_t *);
FLOAT_TYPE Error_ik_k_r( system_t *, parameters_t * );

extern const method_t method_p3m_ik_r;

#endif
