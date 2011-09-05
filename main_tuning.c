#include <stdio.h>

#include "types.h"
#include "common.h"
#include "tuning.h"
#include "error.h"

#include "ewald.h"
#include "p3m-ik.h"
#include "p3m-ik-i.h"
#include "p3m-ad.h"
#include "p3m-ad-i.h"

#include "generate_system.h"

int main( void ) {

  int particles, i;

  FILE *fout = fopen( "timings.dat" , "w" );
  method_t *m[2] = { &method_p3m_ad, &method_p3m_ik };

  system_t *s;

  forces_t *f[2];

  parameters_t *p[2], op;

  data_t *d[2];

  FLOAT_TYPE time[2];

  // strong scaling
  for(particles = 1000; particles < 100000; particles += 1000) {
    printf("Init system with %d particles.", particles);
    #warning no ref forces
    s = generate_system( FORM_FACTOR_RANDOM, particles, 20.0, 1.0 );
    /*
    op.rcut = 0.49*s->length;

    Init_neighborlist( s, &op );
    puts("Calc reference.");
    Calculate_reference_forces( s, &op ); */
   
    for(i=0;i<2;i++) {
      puts("Init.");
      f[i] = Init_forces( s->nparticles );
      puts("Tune");
      p[i] = Tune ( m[i], s, 1e-4 );
      puts("Method init.");
      d[i] = m[i]->Init( s, p[i] );
      puts("Influence function.");
      m[i]->Influence_function( s, p[i], d[i] );
      puts("Init neighborlist");
      Init_neighborlist( s, p[i] );

      puts("Calc...");
      start_timer();
      Calculate_forces ( m[i], s, p[i], d[i], f[i] );
      time [i] = stop_timer();

      puts("Free...");
      Free_data(d[i]);
      Free_forces(f[i]);
      free(p[i]);
    }
    puts("Free systems...");
    Free_system(s);
    fprintf(fout, "%d %e %e\n", particles, time[0], time[1] );
    fflush(fout);
  }
  fclose(fout);
}

