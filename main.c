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

#include "interpol.c"

// Error estimates for p3m

#include "p3m-error.h"

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
    fprintf ( stderr, "usage: %s <positions> <forces> <alpha_min> <alpha_max> <alpha_step> <method>\n", name );
}




int main ( int argc, char **argv ) {
    int methodnr;

    FLOAT_TYPE alphamin,alphamax,alphastep;

    FILE* fout;

    system_t *system;
    method_t method;
    parameters_t parameters;
    data_t *data;
    forces_t *forces;

    error_t error;

    if ( argc != 7 ) {
        usage ( argv[0] );
        return 128;
    }

    // Inits the system and reads particle data and parameters from file.
    system = Daten_einlesen ( &parameters, argv[1] );
    forces = Init_forces(system->nparticles);

    alphamin = atof ( argv[3] );
    alphamax = atof ( argv[4] );
    alphastep = atof ( argv[5] );

    methodnr = atoi ( argv[6] );


    Exakte_Werte_einlesen( system, argv[2] );


    if ( methodnr == method_ewald.method_id )
        method = method_ewald;
#ifdef P3M_IK_H
    else if ( methodnr == method_p3m_ik.method_id )
        method = method_p3m_ik;
#endif
#ifdef P3M_IK_I_H
    else if ( methodnr == method_p3m_ik_i.method_id )
        method = method_p3m_ik_i;
#endif
#ifdef P3M_AD_H
    else if ( methodnr == method_p3m_ad.method_id )
        method = method_p3m_ad;
#endif
#ifdef P3M_AD_I_H
    else if ( methodnr == method_p3m_ad_i.method_id )
        method = method_p3m_ad_i;
#endif
#ifdef GREENS_IK_H
    else if ( methodnr == method_greens_ik.method_id )
        method = method_greens_ik;
#endif
    else {
        fprintf ( stderr, "Method %d not know.", methodnr );
        exit ( 126 );
    }

    if ( ( method.Init == NULL ) || ( method.Influence_function == NULL ) || ( method.Kspace_force == NULL ) ) {
        fprintf ( stderr,"Internal error: Method '%s' (%d) is not properly defined. Aborting.\n", method.method_name, method.method_id );
        exit ( -1 );
    }

    fprintf ( stderr, "Using %s.\n", method.method_name );

    fout = fopen ( "out.dat","w" );

    printf ( "Init" );
    fflush(stdout);
    data = method.Init ( system, &parameters );
    printf ( ".\n" );

    printf ( "Init neighborlist" );
    Init_neighborlist ( system, &parameters, data );
    printf ( ".\n" );

    printf ( "# %8s\t%8s\t%8s\t%8s\t%8s\n", "alpha", "DeltaF", "Estimate", "R-Error", "K-Error" );
    for ( parameters.alpha=alphamin; parameters.alpha<=alphamax; parameters.alpha+=alphastep ) {
        method.Influence_function ( system, &parameters, data );  /* Hockney/Eastwood */

        Calculate_forces ( &method, system, &parameters, data, forces ); /* Hockney/Eastwood */

        error = Calculate_errors ( system, forces );

        if ( method.Error != NULL ) {
            double estimate =  method.Error ( system, &parameters );
            printf ( "%8lf\t%8e\t%8e\t %8e %8e\n", parameters.alpha, error.f / sqrt(system->nparticles) , estimate,
                     error.f_r, error.f_k );
            fprintf ( fout,"% lf\t% e\t% e\n",parameters.alpha,error.f / sqrt(system->nparticles) , estimate );
        } else {
            printf ( "%8lf\t%8e\t na\t%8e\t%8e\n", parameters.alpha,error.f / system->nparticles , error.f_r, error.f_k );
            fprintf ( fout,"% lf\t% e\t na\n",parameters.alpha,error.f / system->nparticles );
        }
#ifdef FORCE_DEBUG
        fprintf ( stderr, "%lf rms %e %e %e\n", parameters.alpha, error.f_v[0], error.f_v[1], error.f_v[2] );
#endif
        fflush ( stdout );
        fflush ( fout );
    }
    fclose ( fout );

    return 0;
}

