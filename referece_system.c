#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "types.h"

#include "p3m-common.h"

// Methods

#include "ewald.h"

#include "interpol.c"

// Utils and IO

#include "io.h"

// Real space part

#include "realpart.h"

// Dipol correction

#include "dipol.h"

// Error calculation

#include "error.h"

// Helper functions for timings

#include "timings.h"

#include "generate_system.h"

// #define WRITE_FORCES

// #define FORCE_DEBUG
// #define CA_DEBUG




void usage ( char *name ) {
    fprintf ( stderr, "usage: %s <nparticles> <box_length> <rcut> <system_file> <forces_file>\n", name );
}




int main ( int argc, char **argv ) {
    int methodnr;

    FLOAT_TYPE boxl;
    int nparticles;

    FILE* fout;

    system_t *system;
    parameters_t parameters;

    forces_t *forces;

    error_t error;

    if ( argc != 6 ) {
        usage ( argv[0] );
        return 128;
    }

    nparticles = atoi(argv[1]);
    boxl = atof(argv[2]);

    parameters.rcut = atof(argv[3]);

    system = generate_system( FORM_FACTOR_RANDOM, nparticles, boxl, 1.0);

    Calculate_reference_forces( system, &parameters );

    Write_system(system, argv[4]);
    Write_exact_forces(system, argv[5]);

    return 0;
}

