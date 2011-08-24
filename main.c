#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "p3m.h"

#include "p3m-common.h"

// Methods


// #include "p3m-ik-i.h"
#include "p3m-ik.h"
// #include "p3m-ad.h"
// #include "p3m-ad-i.h"

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


// #define WRITE_FORCES

// #define FORCE_DEBUG
// #define CA_DEBUG

void Elstat_berechnen ( system_t *, p3m_parameters_t *, const method_t *, p3m_data_t * );

void Elstat_berechnen ( system_t *s, p3m_parameters_t *p, const method_t *m, p3m_data_t *d ) {

    int i, j;

    for ( i=0; i<3; i++ ) {
        memset ( s->f.fields[i]  , 0, s->nparticles*sizeof ( FLOAT_TYPE ) );
        memset ( s->f_k.fields[i], 0, s->nparticles*sizeof ( FLOAT_TYPE ) );
        memset ( s->f_r.fields[i], 0, s->nparticles*sizeof ( FLOAT_TYPE ) );
    }

    start_timer();
    Realteil ( s, p );
    printf ( "Realteil n2: %e\n", stop_timer() );
    start_timer();
    Realpart_neighborlist ( s,p );
    printf ( "Realteil n: %e\n", stop_timer() );
    //  Dipol(s, p);
    start_timer();
    m->Kspace_force ( s, p, d );
    printf ( "k.teil n: %e\n", stop_timer() );
    for ( j=0; j < 3; j++ ) {
#pragma omp parallel for
        for ( i=0; i<s->nparticles; i++ ) {
            s->f.fields[j][i] += s->f_k.fields[j][i] + s->f_r.fields[j][i];
        }
    }
}

void usage ( char *name ) {
    fprintf ( stderr, "usage: %s <positions> <forces> <alpha_min> <alpha_max> <alpha_step> <method>\n", name );
}


void calc_reference_forces ( system_t *s, p3m_parameters_t *p, p3m_data_t *d ) {
  
    int i,j;
    
    p3m_parameters_t op = *p;

    op.alpha = Ewald_compute_optimal_alpha ( s, &op );

    fprintf ( stderr, "Optimal alpha is %lf (error: %e)\n", op.alpha, Ewald_estimate_error ( s, &op ) );

    d = Ewald_init ( s, p );
    Ewald_compute_influence_function ( s, p, d );

    Elstat_berechnen ( s, &op, &method_ewald, d );

    for ( j=0;j<3;j++ )
        for ( i=0;i<s->nparticles;i++ ) {
            s->reference.f.fields[j][i] = s->f.fields[j][i];
            s->reference.f_k.fields[j][i] = s->f_k.fields[j][i];
        }
}


int main ( int argc, char **argv ) {
    int methodnr;

    FLOAT_TYPE alphamin,alphamax,alphastep;

    FILE* fout;

    system_t system;
    method_t method;
    p3m_parameters_t parameters;
    p3m_data_t *data;

    error_t error;

    if ( argc != 7 ) {
        usage ( argv[0] );
        return 128;
    }

    // Inits the system and reads particle data and parameters from file.
    Daten_einlesen ( &system, &parameters, argv[1] );

    // nshift_ausrechnen();

    alphamin = atof ( argv[3] );
    alphamax = atof ( argv[4] );
    alphastep = atof ( argv[5] );

    methodnr = atoi ( argv[6] );

    Interpolationspolynom_berechnen ( parameters.ip ); /* Hockney/Eastwood */

#ifdef WRITE_FORCES
    calc_reference_forces ( argv[2] );
#endif
#ifndef WRITE_FORCES
    Exakte_Werte_einlesen ( &system, argv[2] );
#endif

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
        fprintf ( stderr,"Internal error: Method '%s' (%d) is not correctly defined. Aborting.\n", method.method_name, method.method_id );
        exit ( -1 );
    }

    fprintf ( stderr, "Using %s.\n", method.method_name );

    fout = fopen ( "out.dat","w" );

    printf ( "Init" );
    data = method.Init ( &system, &parameters );
    printf ( ".\n" );

    printf ( "Init neighborlist" );
    Init_neighborlist ( &system, &parameters );
    printf ( ".\n" );

    printf ( "# %8s\t%8s\t%8s\t%8s\t%8s\n", "alpha", "DeltaF", "Estimate", "R-Error", "K-Error" );
    for ( parameters.alpha=alphamin; parameters.alpha<=alphamax; parameters.alpha+=alphastep ) {
        method.Influence_function ( &system, &parameters, data );  /* Hockney/Eastwood */

        Elstat_berechnen ( &system, &parameters, &method, data ); /* Hockney/Eastwood */

        error = Calculate_errors ( &system );

        if ( method.Error != NULL ) {
            double estimate =  method.Error ( &system, &parameters );
            printf ( "%8lf\t%8e\t%8e\t %8e %8e\n", parameters.alpha,error.f, estimate,
                     error.f_r, error.f_k );
            fprintf ( fout,"% lf\t% e\t% e\n",parameters.alpha,error.f, estimate );
        } else {
            printf ( "%8lf\t%8e\t na\t%8e\t%8e\n", parameters.alpha,error.f, error.f_r, error.f_k );
            fprintf ( fout,"% lf\t% e\t na\n",parameters.alpha,error.f );
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

