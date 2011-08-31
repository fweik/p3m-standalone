#include <math.h>
#include <string.h>

#include "error.h"

error_t Calculate_errors(system_t *system, forces_t *f ) {
    int i,j;
    error_t e;
    memset(&e, 0, sizeof(error_t));

    for (i=0; i<system->nparticles; i++) {
        for (j=0;j<3;j++) {
            e.f_v[j] += SQR(system->reference->f->fields[j][i] - f->f->fields[j][i]);

            e.f_k += SQR(system->reference->f_k->fields[j][i] - f->f_k->fields[j][i]);

            e.f_r += SQR(system->reference->f_r->fields[j][i] - f->f_r->fields[j][i]);

            e.f += SQR(f->f->fields[j][i] - system->reference->f->fields[j][i]);
        }
    }
    for (j=0;j<3;j++) {
        e.f_v[j] = sqrt(e.f_v[j]);
    }
    e.f = sqrt(e.f);
    e.f_k = sqrt(e.f_k);
    e.f_r = sqrt(e.f_r);
    return e;
}
