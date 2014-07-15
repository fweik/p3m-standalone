/**    Copyright (C) 2011,2012,2013 Florian Weik <fweik@icp.uni-stuttgart.de>

       This program is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

       This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this program.  If not, see <http://www.gnu.org/licenses/>. **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "types.h"

#include "p3m-common.h"

// Methods

#include "ewald.h"

#include "interpol.h"

// Utils and IO

#include "io.h"

// Real space part

#include "realpart.h"

// Dipol correction

#include "dipol.h"

// Error calculation

#include "error.h"

// Helper functions for timings

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

    system = generate_system( SYSTEM_RANDOM, nparticles, boxl, 1.0);

    Calculate_reference_forces( system, &parameters );

    Write_system(system, argv[4]);
    Write_exact_forces(system, argv[5]);

    return 0;
}

