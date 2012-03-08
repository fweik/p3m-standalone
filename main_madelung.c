#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "types.h"

#include "p3m-common.h"

// Methods


#include "p3m-ik-i.h"
#include "p3m-ik.h"
#include "p3m-ad.h"
#include "p3m-ad-i.h"
#include "greens.h"

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

#include "timings.h"

#include "generate_system.h"

// #define WRITE_FORCES

// #define FORCE_DEBUG
// #define CA_DEBUG


#define MADELUNG -1.747564594633182190636212035544397403485


void usage ( char *name ) {
    fprintf ( stderr, "usage: %s <box_length> <nparticles <rcut>\n", name );
}




int main ( int argc, char **argv ) {
  FLOAT_TYPE ref_prec, rms=0.0;
  int i,j;

    system_t *system;
    method_t method;
    parameters_t parameters;
    data_t *data;
    forces_t *forces;
    FLOAT_TYPE V,M;
    error_t error;

    if ( argc != 4 ) {
        usage ( argv[0] );
        return 128;
    }

    // Inits the system and reads particle data and parameters from file.
    // system = Daten_einlesen ( &parameters, argv[1] );

    system = generate_system( FORM_FACTOR_MADELUNG, atoi(argv[2]), atof(argv[1]), 1.0 );
    parameters.rcut = atof(argv[3]);

    /*    for(i=0; i<system->nparticles;i++) {
      printf("%d pos %e %e %e q %e\n", i, system->p->x[i], system->p->y[i], system->p->z[i], system->q[i]);
      } */

    forces = Init_forces(system->nparticles);

    ref_prec = Calculate_reference_forces( system, &parameters );

    for(i=0;i<system->nparticles;i++) {
      for(j=0;j<3;j++)
	rms += SQR(system->reference->f->fields[j][i]);
    }

    rms = sqrt(rms/system->nparticles);

    V = system->length*system->length*system->length;
    M = 2.0 * system->energy / V;

    printf("rms %e est %e\n", FLOAT_CAST rms, FLOAT_CAST ref_prec );
    printf("energy   %.15f\n", FLOAT_CAST system->energy);
    printf("madelung calc %.40f\n", FLOAT_CAST (2.0 * system->energy / V));
    printf("madelung true %.40f\n", FLOAT_CAST MADELUNG);
    printf("rel. error %e\n", FLOAT_CAST (FLOAT_ABS(MADELUNG - M)/MADELUNG));
}

