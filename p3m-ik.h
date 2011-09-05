#ifndef P3M_IK_H
#define P3M_IK_H

#include "types.h"

void Influence_function_berechnen_ik(system_t*, parameters_t*, data_t*);
void P3M_ik(system_t *, parameters_t *, data_t *, forces_t *);
data_t *Init_ik(system_t*, parameters_t*);
FLOAT_TYPE Error_ik( system_t *, parameters_t *);

extern const method_t method_p3m_ik;

#endif
