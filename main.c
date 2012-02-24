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

void usage ( char *name ) {
    fprintf ( stderr, "usage: %s <positions> <forces> <alpha_min> <alpha_max> <alpha_step> <method>\n", name );
}


static FLOAT_TYPE compute_error_estimate_k(system_t *s, parameters_t *p, FLOAT_TYPE alpha) {
  /* compute the k space part of the error estimate */
  FLOAT_TYPE res, Leni = 1.0/s->length;
  int kmax = p->mesh-1;

  /* Kolafa Perram, eq. 31 */
  res = 2.0 * s->q2 * alpha * Leni * sqrt(1.0/(PI*kmax*s->nparticles))
    * exp(-SQR(PI*kmax/(alpha*s->length)));

  return res;
}


int main ( int argc, char **argv ) {
    int methodnr;

    FLOAT_TYPE alphamin,alphamax,alphastep;

    FILE* fout;

    system_t *system;
    method_t method;
    parameters_t parameters, parameters_ewald;
    data_t *data, *data_ewald;
    forces_t *forces, *forces_ewald;
    char *pos_file = NULL, *force_file = NULL, *out_file = NULL;
    error_t error;
    FLOAT_TYPE length;
    int npart;

    FLOAT_TYPE error_k=0.0, ewald_error_k_est, estimate=0.0, error_k_est;
    int i,j, calc_k_error, calc_est;

    cmd_parameters_t params = { NULL, 0, NULL, 0 };

    add_param( "rcut", ARG_TYPE_FLOAT, ARG_REQUIRED, &(parameters.rcut), &params );
    add_param( "alphamin", ARG_TYPE_FLOAT, ARG_REQUIRED, &alphamin, &params );
    add_param( "alphamax", ARG_TYPE_FLOAT, ARG_REQUIRED, &alphamax, &params );
    add_param( "alphastep", ARG_TYPE_FLOAT, ARG_REQUIRED, &alphastep, &params );
    add_param( "positions", ARG_TYPE_STRING, ARG_OPTIONAL, &pos_file, &params );
    add_param( "forces", ARG_TYPE_STRING, ARG_OPTIONAL, &force_file, &params );
    add_param( "mesh", ARG_TYPE_INT, ARG_REQUIRED, &(parameters.mesh), &params );
    add_param( "cao", ARG_TYPE_INT, ARG_REQUIRED, &(parameters.cao), &params );
    add_param( "method", ARG_TYPE_INT, ARG_REQUIRED, &methodnr, &params );
    add_param( "mc", ARG_TYPE_INT, ARG_OPTIONAL, &P3M_BRILLOUIN, &params );
    add_param( "mc_est", ARG_TYPE_INT, ARG_OPTIONAL, &P3M_BRILLOUIN_TUNING, &params );
    add_param( "error_k", ARG_TYPE_NONE, ARG_OPTIONAL, NULL, &params );
    add_param( "no_estimate", ARG_TYPE_NONE, ARG_OPTIONAL, NULL, &params );
    add_param( "outfile", ARG_TYPE_STRING, ARG_OPTIONAL, &out_file, &params );
    add_param( "particles", ARG_TYPE_INT, ARG_OPTIONAL, &npart, &params );
    add_param( "box", ARG_TYPE_FLOAT, ARG_OPTIONAL, &length, &params );

    parse_parameters( argc - 1, argv + 1, params );

    calc_k_error = param_isset( "error_k", params );
    calc_est = param_isset( "estimate", params );

    parameters.cao3 = parameters.cao*parameters.cao*parameters.cao;
    parameters.ip = parameters.cao - 1;

    parameters_ewald = parameters;

    if( !(param_isset("positions", params) == 1) &&
	!((param_isset("box", params) == 1) && (param_isset("particles", params))) ) {
      puts("Need to provide either 'positions' or 'box' and 'particles'.");
      exit(1);
    }

    if( param_isset("positions", params) == 1) {
    // Inits the system and reads particle data and parameters from file.
      puts("Reading file");
      system = Read_system ( &parameters, pos_file );
      puts("Done.");
    } else {
      puts("Generating system.");
      system = generate_system( FORM_FACTOR_RANDOM, npart, length, 1.0);
      puts("Done.");
    }

    forces = Init_forces(system->nparticles);
    forces_ewald = Init_forces(system->nparticles);

    if(param_isset("forces", params) == 1) {
      printf("Reading reference forces from '%s'.\n", force_file);
      Read_exact_forces( system, force_file );
      puts("Done.");
    } else {
      puts("Calculating reference forces.");
      Calculate_reference_forces( system, &parameters );
      puts("Done.");
    }

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

    if(param_isset("outfile", params) == 1) {
      fout = fopen ( out_file, "w" );      
    } else {
      fout = fopen ( "out.dat","w" );
    }

    printf ( "Init" );
    fflush(stdout);
    data = method.Init ( system, &parameters );
    printf ( ".\n" );

    printf ( "Init Ewald" );
       data_ewald = method_ewald.Init ( system, &parameters_ewald );

    //    printf ( "Init neighborlist" );
    //Init_neighborlist ( system, &parameters, data );
    //printf ( ".\n" );

    printf ( "# %8s\t%8s\t%8s\t%8s\t%8s\n", "alpha", "DeltaF", "Estimate", "R-Error-Est", "K-Error-Est" );
    for ( parameters.alpha=alphamin; parameters.alpha<=alphamax; parameters.alpha+=alphastep ) {
      parameters_ewald.alpha = parameters.alpha;

      method.Influence_function ( system, &parameters, data );  /* Hockney/Eastwood */

      Calculate_forces ( &method, system, &parameters, data, forces ); /* Hockney/Eastwood */

      error_k =0.0;
      if(calc_k_error == 1) {
	for(i=0;i<3;i++) {
	  memset ( forces_ewald->f_k->fields[i], 0, system->nparticles*sizeof ( FLOAT_TYPE ) );
	}
	method_ewald.Influence_function ( system, &parameters_ewald, data_ewald );
	method_ewald.Kspace_force( system, &parameters_ewald, data_ewald, forces_ewald );

	error_k =0.0;
	for (i=0; i<system->nparticles; i++) {
	  for (j=0;j<3;j++) {            
	    error_k   += SQR( forces->f_k->fields[j][i] - forces_ewald->f_k->fields[j][i] );
	  }
	}
	error_k = sqrt(error_k) / sqrt(system->nparticles);
      }

      ewald_error_k_est = compute_error_estimate_k( system, &parameters_ewald, parameters_ewald.alpha);
      error = Calculate_errors ( system, forces );

      if ( method.Error != NULL ) {
	if( calc_est == 0 )
	  estimate = method.Error ( system, &parameters );
	error_k_est = method.Error_k ( system, &parameters);
	printf ( "%8lf\t%8e\t%8e\t %8e %8e\n", FLOAT_CAST parameters.alpha, FLOAT_CAST error.f / sqrt(system->nparticles) , FLOAT_CAST estimate,
		 FLOAT_CAST Realspace_error( system, &parameters ), FLOAT_CAST error_k_est );
	fprintf ( fout,"% lf\t% e\t% e\t% e\t% e\t% e\t% e\n", 
		  FLOAT_CAST parameters.alpha, FLOAT_CAST error.f / sqrt(system->nparticles) , 
		  FLOAT_CAST estimate, Realspace_error( system, &parameters ), 
		  FLOAT_CAST error_k_est, error_k, ewald_error_k_est );
        } else {
            printf ( "%8lf\t%8e\t na\t%8e\t%8e\n", FLOAT_CAST parameters.alpha, FLOAT_CAST error.f / system->nparticles , FLOAT_CAST error.f_r, FLOAT_CAST error.f_k );
            fprintf ( fout,"% lf\t% e\t na\n", FLOAT_CAST parameters.alpha, FLOAT_CAST error.f / system->nparticles );
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

