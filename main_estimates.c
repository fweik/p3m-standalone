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

// Real space part

#include "realpart.h"

// Error calculation

#include "error.h"

// Helper functions for timings

// #define WRITE_FORCES

// #define FORCE_DEBUG
// #define CA_DEBUG




void usage ( char *name ) {
    fprintf ( stderr, "usage: %s <alpha_min> <alpha_max> <alpha_step> <method> <nparticles> <box_length> <cao> <mesh> <rcut> <aliasing_max>\n", name );
}




int main ( int argc, char **argv ) {
    int methodnr;

    FLOAT_TYPE alphamin,alphamax,alphastep, ref_prec;

    FILE* fout;

    system_t system;
    method_t method;
    parameters_t parameters;

    error_t error;

    if ( argc != 11 ) {
        usage ( argv[0] );
        return 128;
    }

    system.nparticles = atoi( argv[5] );
    system.length = atof ( argv[6] );
    system.q2 = system.nparticles;

    parameters.rcut = atof(argv[9]);
    parameters.cao = atoi(argv[7]);
    parameters.ip = parameters.cao - 1;
    parameters.cao3 =  parameters.cao* parameters.cao* parameters.cao;
    parameters.mesh = atoi(argv[8]);

    alphamin = atof ( argv[1] );
    alphamax = atof ( argv[2] );
    alphastep = atof ( argv[3] );

    methodnr = atoi ( argv[4] );

    P3M_BRILLOUIN = atoi(argv[10]);
    P3M_BRILLOUIN_TUNING = P3M_BRILLOUIN;

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
    else {
        fprintf ( stderr, "Method %d not know.", methodnr );
        exit ( 126 );
    }

    if ( ( method.Init == NULL ) || ( method.Influence_function == NULL ) || ( method.Kspace_force == NULL ) ) {
        fprintf ( stderr,"Internal error: Method '%s' (%d) is not properly defined. Aborting.\n", method.method_name, method.method_id );
        exit ( -1 );
    }

    fprintf ( stderr, "Using %s.\n", method.method_name );
    fprintf ( stderr, "System of size %d with box_l %lf.\n", system.nparticles, system.length);

    fout = fopen ( "out.dat","w" );

    printf ( "# %8s\t%8s\n", "alpha", "Estimate" );
    for ( parameters.alpha=alphamin; parameters.alpha<=alphamax; parameters.alpha+=alphastep ) {
        if ( method.Error != NULL ) {
            double estimate =  method.Error ( &system, &parameters );
            printf ( "%8lf\t%8e\n", parameters.alpha, estimate);
            fprintf ( fout,"% lf\t% e\n",parameters.alpha, estimate );
	}
        fflush ( stdout );
        fflush ( fout );
    }
    fclose ( fout );

    return 0;
}

