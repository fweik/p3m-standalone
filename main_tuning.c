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

  FILE *fout = fopen( "timings.dat" , "w" );
  //  FILE *fprec = fopen( "precisions.dat", "w" );
  const method_t *m[4] = { &method_p3m_ik, &method_p3m_ik_i, &method_p3m_ad, &method_p3m_ad_i };

  fprintf(fout, "# particles\t %s\t %s\t %s\t %s\n", m[0]->method_name, m[1]->method_name, m[2]->method_name, m[3]->method_name);

  system_t *s;

  forces_t *f[4];

  parameters_t *p[4], op;

  data_t *d[4];

  //  error_t e[4];

  FLOAT_TYPE time[4];

  // weak scaling

  for(particles = 1000; particles <= 1000; particles += 1000) {
    printf("Init system with %d particles.\n", particles);

    s = generate_system( FORM_FACTOR_RANDOM, particles, my_power(particles / 10000, 1.0/3.0) * 20.0, 1.0 );
 
    //op.rcut = 0.49*s->length;

    //    puts("Calc reference.");
    //Calculate_reference_forces( s, &op );

    // methods are independent of each other

    for(i=0;i<4;i++) {
      puts("Init.");
      f[i] = Init_forces( s->nparticles );
      puts("Tune");
      p[i] = Tune ( m[i], s, 1e-4 );
      
      puts("Method init.");
      d[i] = m[i]->Init( s, p[i] );
      puts("Influence function.");
      m[i]->Influence_function( s, p[i], d[i] );
      puts("Init neighborlist");
      Init_neighborlist( s, p[i], d[i] );

      puts("Calc...");

      time[i] = MPI_Wtime();      
      for(j=0;j<50;j++)
	Calculate_forces ( m[i], s, p[i], d[i], f[i] );
      time[i] = MPI_Wtime() - time[i];

      //e[i] = Calculate_errors( s, f[i] );

      Free_neighborlist(d[i]);

      puts("Free...");
      Free_data(d[i]);
      Free_forces(f[i]);
      fftw_free(p[i]);
    }

    puts("Free systems...");
    Free_system(s);
    fprintf(fout, "%d %e %e %e %e\n", particles, time[0], time[1], time[2], time[3] );
    //    fprintf(fprec, "%d %e %e %e %e\n", particles, e[0].f / particles, e[1].f / particles, e[2].f / particles, e[2].f / particles );
    fflush(fout);
    //fflush(fprec);
  }
  fclose(fout);
  //  fclose(fprec);

  return 0;
}

