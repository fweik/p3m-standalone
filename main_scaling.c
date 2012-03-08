#include <stdio.h>
#include <mpi.h>
#include <math.h>

#include "types.h"
#include "common.h"
#include "tuning.h"
#include "error.h"
#include "p3m-common.h"

#include "ewald.h"
#include "p3m-ik.h"
#include "p3m-ik-i.h"
#include "p3m-ad.h"
#include "p3m-ad-i.h"

#include "realpart.h"

#include "generate_system.h"

int main( void ) {

  int particles, i, j;

  const method_t *m[4] = { &method_p3m_ik, &method_p3m_ik_i, &method_p3m_ad, &method_p3m_ad_i };

  system_t *s;

  forces_t *f[4];

  parameters_t *p[4];

  data_t *d[4];

  FLOAT_TYPE time[4];

  FILE *fout = fopen("scaling.d.dat", "w");

  int incr = 100;

  //  FUKW *ferror = fropen("error.d.dat", "w");

  // weak scaling

  for(particles = 100; particles <= 10000; particles += incr) {
    fprintf(stderr, "Init system with %d particles", particles);

    if(particles >= 1000)
      incr = 1000;
    if(particles >= 10000)
      incr = 10000;
    if(particles >= 100000)
      incr = 100000;

    s = generate_system( FORM_FACTOR_RANDOM, particles, pow (particles / 1000.0, 1.0/3.0) * 10.0, 1.0 );
 
    fprintf(stderr, ".\n");


     
#pragma omp barrier
    fprintf( stderr, "starting...\n" );

    fprintf(fout, "%d ", particles);


    // methods are independent of each other
#pragma omp parallel for private(j)
    for(i=0;i<4;i++) {
      int mesh;

      f[i] = Init_forces( s->nparticles );

      p[i] = Init_array( 1, sizeof ( parameters_t ) );

      fprintf( stderr, "Method '%s'\n", m[i]->method_name);

      mesh = (int) pow( particles, 1.0 / 3.0 );

      mesh = ( mesh >= 8 ) ? mesh : 8;

      p[i]->mesh = mesh;
      p[i]->alpha = 1.0;
      p[i]->cao = 7;
      p[i]->cao3 = 7*7*7;
      p[i]->ip = 6;
      p[i]->rcut = 3.0;

      d[i] = m[i]->Init( s, p[i] );

      m[i]->Influence_function( s, p[i], d[i] );

      time[i] = MPI_Wtime();      
      m[i]-> Kspace_force( s, p[i], d[i], f[i] );
      fprintf( fout, " %e", MPI_Wtime() - time[i]);
      
      Free_data(d[i]);
      Free_forces(f[i]);
      FFTW_FREE(p[i]);

      fprintf(stderr, "done...\n");

      fflush(stderr);
      fflush(fout);
      //  fflush(ferror);
    }
    fprintf(fout, "\n");
    #pragma omp barrier

    Free_system(s);
  }
  fclose(fout);
  return 0;
}

