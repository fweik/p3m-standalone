#ifndef GREENS_IK_H
#define GREENS_IK_H

#include "types.h"

void Greens_function(system_t*, parameters_t*, data_t*);
void Greens_kspace_ik(system_t *, parameters_t *, data_t *, forces_t *);
data_t *Init_greens_ik(system_t*, parameters_t*);

extern const method_t method_greens_ik;

#endif
