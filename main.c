#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "types.h"

#include "p3m-common.h"

// Methods

#include "p3m-ik-i.h"
#include "p3m-ik.h"
#include "p3m-ad.h"
#include "p3m-ad-i.h"

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

#include "statistics.h"
#include "common.h"

#include "inhomogenous_error.h"

// #define WRITE_FORCES

// #define FORCE_DEBUG
// #define CA_DEBUG

static FLOAT_TYPE compute_error_estimate_k(system_t *s, parameters_t *p, FLOAT_TYPE alpha) {
  /* compute the k space part of the error estimate */
  FLOAT_TYPE res, Leni = 1.0/s->length;
  int kmax = p->mesh-1;

  /* Kolafa Perram, eq. 31 */
  res = 2.0 * s->q2 * alpha * Leni * SQRT(1.0/(PI*kmax*s->nparticles))
    * EXP(-SQR(PI*kmax/(alpha*s->length)));

  return res;
}


int main ( int argc, char **argv ) {
    int methodnr;

    FLOAT_TYPE alphamin,alphamax,alphastep;
    FLOAT_TYPE alpha;
    FLOAT_TYPE wtime;

    FILE* fout;
    
    system_t *system;
    method_t method;
    parameters_t parameters, parameters_ewald;
    data_t *data, *data_ewald;
    forces_t *forces, *forces_ewald;
    char *pos_file = NULL, *force_file = NULL, *out_file = NULL, *ref_out = NULL, *sys_out = NULL, *rdf_file = NULL, *vtf_file = NULL, *cdf_file = NULL;
    error_t error;
    FLOAT_TYPE length, prec;
    int npart;
    FLOAT_TYPE charge;
    int form_factor;
    FLOAT_TYPE rdf_min, rdf_max;
    int rdf_bins;

    int inhomo_error_mesh = 64;
    int inhomo_error_cao = 5;

    FLOAT_TYPE error_k=0.0, ewald_error_k_est, estimate=0.0, error_k_est;
    int i,j, calc_k_error, calc_est;

    int error_map_mesh=64, error_map_cao=1;

    #ifdef _OPENMP
    int nthreads;
    #endif

    cmd_parameters_t params = { NULL, 0, NULL, 0 };

    add_param( "rcut", ARG_TYPE_FLOAT, ARG_REQUIRED, &(parameters.rcut), &params );
    add_param( "alphamin", ARG_TYPE_FLOAT, ARG_OPTIONAL, &alphamin, &params );
    add_param( "alphamax", ARG_TYPE_FLOAT, ARG_OPTIONAL, &alphamax, &params );
    add_param( "alphastep", ARG_TYPE_FLOAT, ARG_OPTIONAL, &alphastep, &params );
    add_param( "alpha", ARG_TYPE_FLOAT, ARG_OPTIONAL, &alpha, &params );
    add_param( "positions", ARG_TYPE_STRING, ARG_OPTIONAL, &pos_file, &params );
    add_param( "error_k", ARG_TYPE_NONE, ARG_OPTIONAL, NULL, &params );
    add_param( "forces", ARG_TYPE_STRING, ARG_OPTIONAL, &force_file, &params );
    add_param( "mesh", ARG_TYPE_INT, ARG_REQUIRED, &(parameters.mesh), &params );
    add_param( "cao", ARG_TYPE_INT, ARG_REQUIRED, &(parameters.cao), &params );
    add_param( "method", ARG_TYPE_INT, ARG_REQUIRED, &methodnr, &params );
    add_param( "mc", ARG_TYPE_INT, ARG_OPTIONAL, &P3M_BRILLOUIN, &params );
    add_param( "mc_est", ARG_TYPE_INT, ARG_OPTIONAL, &P3M_BRILLOUIN_TUNING, &params );
    add_param( "no_estimate", ARG_TYPE_NONE, ARG_OPTIONAL, NULL, &params );
    add_param( "outfile", ARG_TYPE_STRING, ARG_OPTIONAL, &out_file, &params );
    add_param( "particles", ARG_TYPE_INT, ARG_OPTIONAL, &npart, &params );
    add_param( "box", ARG_TYPE_FLOAT, ARG_OPTIONAL, &length, &params );
    add_param( "tune", ARG_TYPE_NONE, ARG_OPTIONAL, NULL, &params );
    add_param( "prec", ARG_TYPE_FLOAT, ARG_OPTIONAL, &prec, &params );
    add_param( "reference_out", ARG_TYPE_STRING, ARG_OPTIONAL, &ref_out, &params );
    add_param( "system_out", ARG_TYPE_STRING, ARG_OPTIONAL, &sys_out, &params );
    add_param( "verlet_lists", ARG_TYPE_NONE, ARG_OPTIONAL, NULL, &params );
    add_param( "charge", ARG_TYPE_FLOAT, ARG_OPTIONAL, &charge, &params );
    add_param( "system_type", ARG_TYPE_INT, ARG_OPTIONAL, &form_factor, &params );
    add_param( "rdf", ARG_TYPE_STRING, ARG_OPTIONAL, &rdf_file, &params );
    add_param( "rdf_bins", ARG_TYPE_INT, ARG_OPTIONAL, &rdf_bins, &params );
    add_param( "rdf_rmin", ARG_TYPE_FLOAT, ARG_OPTIONAL, &rdf_min, &params );
    add_param( "rdf_rmax", ARG_TYPE_FLOAT, ARG_OPTIONAL, &rdf_max, &params );
    add_param( "cdf", ARG_TYPE_STRING, ARG_OPTIONAL, &cdf_file, &params );
    add_param( "no_calculation", ARG_TYPE_NONE, ARG_OPTIONAL, NULL, &params );
    add_param( "vtf_file", ARG_TYPE_STRING, ARG_OPTIONAL, &vtf_file, &params );
    add_param( "rdf_species", ARG_TYPE_NONE, ARG_OPTIONAL, NULL, &params );
    add_param( "no_reference_force", ARG_TYPE_NONE, ARG_OPTIONAL, NULL, &params);
    #ifdef _OPENMP
    add_param( "threads", ARG_TYPE_INT, ARG_OPTIONAL, &nthreads, &params );
    #endif
    add_param( "inhomo_mesh", ARG_TYPE_INT, ARG_OPTIONAL, &inhomo_error_mesh, &params);
    add_param( "inhomo_cao", ARG_TYPE_INT, ARG_OPTIONAL, &inhomo_error_cao, &params);
    add_param( "error_map", ARG_TYPE_NONE, ARG_OPTIONAL, NULL, &params);
    add_param( "error_map_mesh", ARG_TYPE_INT, ARG_OPTIONAL, &error_map_mesh, &params);
    add_param( "error_map_cao", ARG_TYPE_INT, ARG_OPTIONAL, &error_map_cao, &params);

    parse_parameters( argc - 1, argv + 1, params );

    calc_k_error = param_isset( "error_k", params );
    calc_est = param_isset( "estimate", params );

    parameters.cao3 = parameters.cao*parameters.cao*parameters.cao;
    parameters.ip = parameters.cao - 1;
    parameters.alpha = 0.0;
    parameters_ewald = parameters;

#ifdef _OPENMP
    if(param_isset("threads", params)) {
      omp_set_num_threads(nthreads);
    }
    printf("OpenMP: Using up to %d threads.\n", omp_get_max_threads( ));
#endif

    if(!(param_isset("alphamin", params) && param_isset("alphamax", params) && param_isset("alphastep", params)) && !param_isset("alpha", params)) {
      puts("Need to provide either alpha-range (alphamin, alphamax, alphastep) or alpha.");
      exit(1);
    }

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
      if( !(param_isset("charge", params) == 1)) {
	charge=1.0;
      }
      if(param_isset("system_type", params)) {
	printf("Using system type %d\n", form_factor);
	system = generate_system( form_factor, npart, length, charge);
      } else {
	system = generate_system( FORM_FACTOR_RANDOM, npart, length, charge);
      }
      puts("Done.");
    }

    /* inhomo_error(system, NULL, 100); */

    if( param_isset("vtf_file", params) == 1) 
      write_vtf( vtf_file, system );

    if( param_isset("rdf", params) == 1) {
      puts("Calculating RDF");
      if( param_isset("rdf_bins", params) == 0)
	rdf_bins = 100;
      if( param_isset("rdf_rmin", params) == 0)
	rdf_min = 0.0;
      if( param_isset("rdf_rmax", params) == 0)
	rdf_max = system->length/2;
      printf("Using %d bins, %lf <= r <= %lf\n", rdf_bins, rdf_min, rdf_max);
      int bins = rdf_bins;
      FLOAT_TYPE *rdf = radial_distribution(rdf_min, rdf_max, rdf_bins, system);
      FLOAT_TYPE *c;
      /* FLOAT_TYPE *rdf_sym = Init_array( 2*bins-1, 2*sizeof(FLOAT_TYPE)); */
      /* c = low_pass_forward( bins, rdf, 0.3); */
      /* c = low_pass_backward(bins, c, 0.3); */
      /* rdf = c; */
      /* for(int i = 0; i < 2*bins; i++) { */
      /* 	rdf_sym[i] = c[i]; */
      /* } */
      /* for(int i = bins; i < 2*bins-1; i++) { */
      /* 	rdf_sym[2*i] =  c[2*bins - 2] + (i-bins)*(c[2] - c[0]); */
      /* 	rdf_sym[2*i+1] = c[4*bins - 2*i - 1]; */
      /* } */


      FILE *rdf_out = fopen(rdf_file, "w");
      /* FILE *c_out = fopen("c_fft.dat", "w"); */

      /* rshif_array(2*(2*bins-1), c, 2*bins); */

      for(int i = 0; i<bins; i++)
	fprintf(rdf_out, "%e %e\n", FLOAT_CAST rdf[2*i], FLOAT_CAST rdf[2*i+1]);

      /* for(int i = 0; i<2*bins-1; i++) */
      /* 	/\* fprintf(c_out, "%e %e\n", FLOAT_CAST rdf_sym[2*i], FLOAT_CAST rdf_sym[2*i+1] );  *\/ */
      /* 	fprintf(c_out, "%d %e %e\n", i, FLOAT_CAST c[2*i], FLOAT_CAST c[2*i+1] ); */

      /* fclose(c_out); */
      fclose(rdf_out);
      fftw_free(rdf);
      /* fftw_free(c); */
      puts("Done.");
    }

    if( param_isset("rdf_species", params) == 1) {
      radial_distribution_species(0.0, 3.0, 200, system);
    }

    forces = Init_forces(system->nparticles);
    forces_ewald = Init_forces(system->nparticles);

    if(param_isset("reference_out", params)) 
      {
	printf("Minimal distance: %.*f\n", DIGITS, FLOAT_CAST Min_distance( system ));
	puts("Calculating reference forces.");
	printf("Reference precision %e\n.", FLOAT_CAST Calculate_reference_forces( system, &parameters ));
	puts("Done.");
	printf("Writing reference forces to '%s'\n", ref_out);
	Write_exact_forces(system, ref_out);
	puts("Done.");
      }

    if(param_isset("system_out", params)) 
      {
	printf("Writing system to '%s'\n", sys_out);
	Write_system(system, sys_out);
	puts("Done.");
      }

    if(param_isset("no_calculation", params) == 1)
      return 0;

    if(param_isset("forces", params) == 1) {
      printf("Reading reference forces from '%s'.\n", force_file);
      Read_exact_forces( system, force_file );
      puts("Done.");
    } else {
      if(param_isset("no_reference_force", params) !=1) {
	puts("Calculating reference forces.");
	printf("Reference precision %e\n.", FLOAT_CAST Calculate_reference_forces( system, &parameters ));
	puts("Done.");
      } else {
	puts("Skipping reference force calculation.");
      }
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
    else if ( methodnr == method_p3m_ad_i.method_id ) {
        method = method_p3m_ad_i;
    }
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
    printf ( ".\n" );

    /* printf ( "Init neighborlist" ); */
    /* Init_neighborlist ( system, &parameters, data ); */
    /* printf ( ".\n" ); */

    FLOAT_TYPE gen_err_dip, gen_err;

    printf ( "# %8s\t%8s\t%8s\t%8s\t%8s\n", "alpha", "DeltaF", "Estimate", "R-Error-Est", "K-Error-Est Generic-K-Space-err" );
    for ( parameters.alpha=alphamin; parameters.alpha<=alphamax; parameters.alpha+=alphastep ) {
      parameters_ewald.alpha = parameters.alpha;

      method.Influence_function ( system, &parameters, data );  /* Hockney/Eastwood */

      wtime = MPI_Wtime();

      Calculate_forces ( &method, system, &parameters, data, forces ); /* Hockney/Eastwood */

      wtime = MPI_Wtime() - wtime;

      error_k =0.0;
      if(calc_k_error == 1) {
	puts("kerror");
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
	error_k = SQRT(error_k) / SQRT(system->nparticles);
	if(param_isset("error_map", params)) {
	  puts("Writing error map.");
	  FLOAT_TYPE *error_map;
	  error_map = Error_map(system, forces, forces_ewald, error_map_mesh, error_map_cao);
	  FILE *error_out = fopen("error_map.dat", "w");
	  int nx,ny,nz;
	  int ind;
	  for (nx=0; nx<error_map_mesh; nx++) {
	    for (ny=0; ny<error_map_mesh; ny++) {
	      for (nz=0; nz<error_map_mesh; nz++) {
		ind = 2*((error_map_mesh*error_map_mesh*nx) + error_map_mesh*(ny) + (nz));
		fprintf(error_out, "%d %d %d %e\n", nx, ny, nz, FLOAT_CAST error_map[ind]);
	      }
	    }
	  }
	  FFTW_FREE(error_map);
	}
      }

      ewald_error_k_est = compute_error_estimate_k( system, &parameters_ewald, parameters_ewald.alpha);
      error = Calculate_errors ( system, forces );

      if ( method.Error != NULL ) {
	if( calc_est == 0 )
	  estimate = method.Error ( system, &parameters );
	error_k_est = method.Error_k ( system, &parameters);

	FLOAT_TYPE err_inhomo = Generic_error_estimate_inhomo(system, &parameters, inhomo_error_mesh, inhomo_error_cao, P3M_BRILLOUIN);
	FLOAT_TYPE rs_error = Realspace_error( system, &parameters );

	/* printf("Q_uncorr %e, Q_corr %e, Q_nonfluc %e\n", Q_uncorr, Q_corr, Q_nonfluc); */

	printf ( "%8lf\t%8e\t%8e\t %8e %8e\t %8e sec\t %8e\n", FLOAT_CAST parameters.alpha, FLOAT_CAST (error.f / SQRT(system->nparticles)) , FLOAT_CAST estimate,
		 FLOAT_CAST rs_error , FLOAT_CAST error_k_est, FLOAT_CAST wtime, FLOAT_CAST err_inhomo );

	/* printf ( "%8lf\t%8e\t%8e\t %8e %8e\t %8e sec\n", FLOAT_CAST parameters.alpha, FLOAT_CAST (error.f / SQRT(system->nparticles)) , FLOAT_CAST estimate, */
	/* 	 FLOAT_CAST rs_error , FLOAT_CAST error_k_est, FLOAT_CAST wtime ); */

	fprintf ( fout,"% lf\t% e\t% e\t% e\t% e\t% e\t% e\t% e\n", 
		  FLOAT_CAST parameters.alpha, FLOAT_CAST (error.f / SQRT(system->nparticles)) , 
		  FLOAT_CAST estimate, FLOAT_CAST Realspace_error( system, &parameters ), 
		  FLOAT_CAST error_k_est, FLOAT_CAST error_k, FLOAT_CAST ewald_error_k_est, FLOAT_CAST err_inhomo);
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

