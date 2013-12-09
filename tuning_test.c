#include <string.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "types.h"
#include "generate_system.h"
#include "tuning.h"
#include "error.h"
#include "p3m-common.h"

#include "p3m-ik.h"
#include "p3m-ik-i.h"
#include "p3m-ad.h"
#include "p3m-ad-i.h"
#include "p3m-ik-real.h"
#include "p3m-ad-real.h"

int main(int argc, char **argv) {
  system_t *s;
  parameters_t p;
  int start, stop, step;  
  FLOAT_TYPE prec;
  method_t methods[6] = { method_p3m_ik, method_p3m_ik_r, method_p3m_ik_i, method_p3m_ad, method_p3m_ad_i, method_p3m_ad_r };
  FILE *f[6];
  FILE *sys_params;
  FLOAT_TYPE box;
  FLOAT_TYPE rcut;
  char filename_buffer[256];
  FLOAT_TYPE density;
  FLOAT_TYPE charge;
  
  start = atoi(argv[1]);
  stop = atoi(argv[2]);
  step = atoi(argv[3]);
  prec = atof(argv[4]);
  rcut = atof(argv[5]);
  density = atof(argv[6]);
  charge = atof(argv[7]);

  for(int i = 0; i < 6; i++) {
    sprintf(filename_buffer, "%s-%d-to-%d-prec-%.1e-rcut-%1.1lf-dens-%1.1lf-q-%1.1lf.dat", methods[i].method_name_short,
	    start, stop, prec, rcut, density, charge);
    f[i] = fopen( filename_buffer, "w");
    fprintf(f[i], "# %s\n", methods[i].method_name);
    fprintf(f[i], "# number_of_part mesh cao alpha time\n");
  }

  sprintf(filename_buffer, "%s-%d-to-%d-prec-%.1e-rcut-%1.1lf-dens-%1.1lf-q-%1.1lf.dat", "params",
	  start, stop, prec, rcut, density, charge);
  sys_params = fopen( filename_buffer, "w");

  fprintf(sys_params, "# npart density box\n");

  for(int i = start; i <= stop; i+= step) {
    printf("Tuning for %d particles.\n", i);
    box = pow((double)(i)/density, 0.3333);
    printf("density %lf (%lf), box %lf npart %d \n", FLOAT_CAST (double)(i)/(box*box*box), FLOAT_CAST density, FLOAT_CAST box, i);    
    fprintf(sys_params, "%d %lf %lf\n", i, FLOAT_CAST (double)(i)/(box*box*box), FLOAT_CAST box);
    fflush(sys_params);

    s = generate_system( SYSTEM_RANDOM, i, box, charge);

    for(int j = 0; j < 6; j++) {
      memset(&p, 0, sizeof(parameters_t));
      p.rcut = rcut;

      timing_t timing;

      timing = Tune( methods+j, s, &p, prec);
      if( timing.avg < 0.0) {
	puts("Tuning failed.");
	continue;
      }
      printf("\t%s:\n", methods[j].method_name);
      printf("\t\tmesh %d cao %d time %lf +/- %lf\n", p.mesh, p.cao, timing.avg, timing.sgm);

      forces_t *forces = Init_forces(i);
  
      data_t *d = methods[j].Init( s, &p);
  
      Free_data(d);
      Free_forces(forces);

      

      fprintf(f[j], "%d %d %d %lf %e %e\n", i, p.mesh, p.cao, p.alpha, timing.avg, timing.sgm);
      fflush(f[j]);
    }

    Free_system(s);
  }

  fclose(sys_params);
  for(int i = 0; i < 6; i++)
    fclose(f[i]);
}
