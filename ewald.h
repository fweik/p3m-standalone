#ifndef EWALD_H
#define EWALD_H

#include "types.h"

#ifndef r_ind
#define r_ind(A,B,C) ((A)*d->mesh*d->mesh + (B)*d->mesh + (C))
#endif

data_t *Ewald_init(system_t *, parameters_t *);
void Ewald_k_space(system_t *, parameters_t *, data_t *, forces_t *);
void Ewald_compute_influence_function(system_t *, parameters_t *, data_t *);
FLOAT_TYPE Ewald_compute_optimal_alpha(system_t *, parameters_t *);
FLOAT_TYPE Ewald_estimate_error(system_t *, parameters_t *);

extern const method_t method_ewald;

#endif
