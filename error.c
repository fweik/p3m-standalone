#include <math.h>
#include <string.h>

#include "error.h"

error_t Calculate_errors(system_t *system, forces_t *f ) {
    int i,j;
    error_t e;
    memset(&e, 0, sizeof(error_t));

    for (i=0; i<system->nparticles; i++) {
        for (j=0;j<3;j++) {
            e.f   += SQR(  f->f->fields[j][i] - system->reference->f->fields[j][i] );
        }
    }

    e.f = sqrt(e.f);
    return e;
}
