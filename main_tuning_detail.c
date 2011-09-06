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

#ifdef DETAILED_TIMINGS
int __detailed_timings=0;
#endif

int main( void ) {

  int particles, i, j;

  const method_t *m[4] = { &method_p3m_ik, &method_p3m_ik_i, &method_p3m_ad, &method_p3m_ad_i };

  system_t *s;

  forces_t *f[4];

  parameters_t *p[4], op;

  data_t *d[4];

  error_t e[4];

  FLOAT_TYPE time[4];

  // weak scaling

  for(particles = 10000; particles <= 15000; particles += 1000) {

    s = generate_system( FORM_FACTOR_RANDOM, particles, pow(particles / 1000, 1.0/3.0) * 20.0, 1.0 );
 
    op.rcut = 0.49*s->length;

    Calculate_reference_forces( s, &op );

#pragma omp barrier
    fprintf( stderr, "starting...\n" );
    printf("%d", particles);

    // methods are independent of each other
#pragma omp parallel for private(j)
    for(i=0;i<4;i++) {

      f[i] = Init_forces( s->nparticles );

      p[i] = Tune ( m[i], s, 1e-4 );
      
      fprintf(stderr, "tuned...\n");

      d[i] = m[i]->Init( s, p[i] );

      m[i]->Influence_function( s, p[i], d[i] );

      Init_neighborlist( s, p[i], d[i] );

      __detailed_timings = 1;

      time[i] = MPI_Wtime();      
	Calculate_forces ( m[i], s, p[i], d[i], f[i] );

      printf(" %e", MPI_Wtime() - time[i]);

      __detailed_timings = 0;

      e[i] = Calculate_errors( s, f[i] );

      Free_neighborlist(d[i]);

      Free_data(d[i]);
      Free_forces(f[i]);
      fftw_free(p[i]);

      fprintf(stderr, "done...\n");
      fflush(stdin);
    }
    printf("\n");
    #pragma omp barrier

    Free_system(s);
  }
  return 0;
}

