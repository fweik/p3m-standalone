#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <string.h>

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

int main( int argc, char *argv[] ) {

  int particles, i, j;
  char filename[] = "det.dat";

  const FLOAT_TYPE rcut = 3.0;

  const method_t *m[4] = { &method_p3m_ik, &method_p3m_ad, &method_p3m_ik_i, &method_p3m_ad_i };

  FILE *fout = fopen( "timings.dat" , "w" );
  FILE *fmesh = fopen( "mesh.dat", "w" );
  FILE *falpha = fopen( "alpha.dat", "w" );
  FILE *fcao = fopen( "cao.dat", "w" );
  FILE *fprec = fopen( "prec.dat", "w" );
#ifdef __detailed_timings
  FILE *timings[4];
  for(i=0;i<4;i++) {
    filename[0] = 48 + m[i]->method_id;
    timings[i] = fopen( filename, "w");
  }
#endif

  fprintf(fout, "# rcut %e\n", rcut);
  fprintf(fout, "# particles\t %s\t %s\t %s\t %s\n", m[0]->method_name, m[1]->method_name, m[2]->method_name, m[3]->method_name);
  fprintf(fmesh, "# rcut %e\n", rcut);
  fprintf(fmesh, "# particles\t %s\t %s\t %s\t %s\n", m[0]->method_name, m[1]->method_name, m[2]->method_name, m[3]->method_name);
  fprintf(falpha, "# rcut %e\n", rcut);
  fprintf(falpha, "# particles\t %s\t %s\t %s\t %s\n", m[0]->method_name, m[1]->method_name, m[2]->method_name, m[3]->method_name);
  fprintf(fcao, "# rcut %e\n", rcut);
  fprintf(fcao, "# particles\t %s\t %s\t %s\t %s\n", m[0]->method_name, m[1]->method_name, m[2]->method_name, m[3]->method_name);
  fprintf(fprec, "# rcut %e\n", rcut);
  fprintf(fprec, "# particles\t %s\t %s\t %s\t %s\n", m[0]->method_name, m[1]->method_name, m[2]->method_name, m[3]->method_name);

  system_t *s;

  forces_t *f[4];

  parameters_t *p[4], op;

  data_t *d[4];

  //  error_t e[4];

  FLOAT_TYPE time[4];
  FLOAT_TYPE precision=1e-4;
  double boxl;

  int part_min, part_max, part_step;

  // weak scaling

  for(i=0;i<4;i++) {
    p[i] = Init_array( 1, sizeof(parameters_t) );
    memset( p[i], 0, sizeof(parameters_t) );
    p[i]->rcut = rcut;
  }

  if(argc == 5) {
    part_min = atoi(argv[1]);
    part_max = atoi(argv[2]);
    part_step = atoi(argv[3]);
    precision = atof(argv[4]);
  } else {
    part_min = 1000;
    part_max = 10000;
    part_step = 1000;
  }

  printf("Timing for %d to %d particles in increments of %d.\n", part_min, part_max, part_step);

  for(particles = part_min; particles <= part_max; particles += part_step) {
    boxl = pow(particles / 100.0, 1.0/3.0)*20.0;
    printf("Init system with %d particles.\n", particles);
    printf("box %e\n", boxl);

    s = generate_system( FORM_FACTOR_RANDOM, particles, boxl, 1.0 );
 
    //    op.rcut = 0.49*s->length;

    //        puts("Calc reference.");
    //Calculate_reference_forces( s, &op );

    // methods are independent of each other

    for(i=0;i<4;i++) {
      puts("Init.");
      f[i] = Init_forces( s->nparticles );
      puts("Tune");

      P3M_BRILLOUIN = 0;
      P3M_BRILLOUIN_TUNING = 0;
      if(m[i]->flags & METHOD_FLAG_ad)
	P3M_BRILLOUIN_TUNING = 1;

      if(m[i]->flags & METHOD_FLAG_interlaced) {
	puts("Interlaced method.");
	P3M_BRILLOUIN = 1;
	P3M_BRILLOUIN_TUNING = 1;
      }

      memset( p[i], 0, sizeof(parameters_t) );
      p[i]->rcut = 3.0;

      if( Tune ( m[i], s, p[i], precision ) != 0)
	continue;
      
      puts("Method init.");
      d[i] = m[i]->Init( s, p[i] );
      puts("Influence function.");
      m[i]->Influence_function( s, p[i], d[i] );

      puts("Calc...");

      time[i] = MPI_Wtime();      
      m[i]->Kspace_force( s, p[i], d[i], f[i] );
      time[i] = MPI_Wtime() - time[i];

#ifdef __detailed_timings
      fprintf(timings[i], "%d %e %e %e %e\n", particles, t_charge_assignment[i], t_force_assignment[i], t_convolution[i], t_fft[i]);
      fflush(timings[i]);	      
#endif


      puts("Free...");
      Free_data(d[i]);
      Free_forces(f[i]);
    }

    puts("Free systems...");
    Free_system(s);
    fprintf(fout, "%d %e %e %e %e\n", particles, time[0], time[1], time[2], time[3] );
    fprintf(fmesh, "%d %d %d %d %d\n", particles, p[0]->mesh, p[1]->mesh, p[2]->mesh, p[3]->mesh );
    fprintf(falpha, "%d %e %e %e %e\n", particles, p[0]->alpha, p[1]->alpha, p[2]->alpha, p[3]->alpha );
    fprintf(fcao, "%d %d %d %d %d\n", particles, p[0]->cao, p[1]->cao, p[2]->cao, p[3]->cao );

    //    fprintf(fprec, "%d %e %e %e %e\n", particles, e[0].f / particles, e[1].f / particles, e[2].f / particles, e[2].f / particles );
    fflush(fout);
    fflush(fmesh);
    fflush(falpha);
    fflush(fcao);

    //fflush(fprec);
  }
  fclose(fout);
  //  fclose(fprec);

  return 0;
}

