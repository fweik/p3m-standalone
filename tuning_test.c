#include <string.h>
#include <stdio.h>
#include <math.h>

#include "types.h"
#include "generate_system.h"
#include "tuning.h"
#include "error.h"

#include "p3m-ik.h"
#include "p3m-ik-i.h"
#include "p3m-ad.h"
#include "p3m-ad-i.h"
#include "p3m-ik-real.h"

int main(int argc, char **argv) {
  system_t *s;
  parameters_t p;
  int start, stop, step;  
  FLOAT_TYPE prec;
  method_t methods[5] = { method_p3m_ik, method_p3m_ik_i, method_p3m_ad, method_p3m_ad_i, method_p3m_ik_r };
  FILE *f[5];
  FLOAT_TYPE time;
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

  for(int i = 0; i < 5; i++) {
    sprintf(filename_buffer, "%s-%d-to-%d-prec-%.1e-rcut-%1.1lf-dens-%1.1lf-q-%1.1lf.dat", methods[i].method_name_short,
	    start, stop, prec, rcut, density, charge);
    f[i] = fopen( filename_buffer, "w");
    fprintf(f[i], "# %s\n", methods[i].method_name);
    fprintf(f[i], "# number_of_part mesh cao alpha time\n");
  }
  for(int i = start; i <= stop; i+= step) {
    printf("Tuning for %d particles.\n", i);
    box = pow((double)(i)/density, 0.3333);
    printf("density %lf (%lf), box %lf npart %d \n", FLOAT_CAST (double)(i)/(box*box*box), FLOAT_CAST density, FLOAT_CAST box, i);    
    s = generate_system( FORM_FACTOR_RANDOM, i, box, charge);

    for(int j = 0; j < 5; j++) {
      memset(&p, 0, sizeof(parameters_t));
      p.rcut = rcut;

      time = Tune( methods+j, s, &p, prec);
      if( time < 0.0) {
	puts("Tuning failed.");
	continue;
      }
      printf("\t%s:\n", methods[j].method_name);
      printf("\t\tmesh %d cao %d time %lf\n", p.mesh, p.cao, time);

      fprintf(f[j], "%d %d %d %lf %e\n", i, p.mesh, p.cao, p.alpha, time);
      fflush(f[j]);
    }

    Free_system(s);
  }

  for(int i = 0; i < 5; i++)
    fclose(f[i]);
}
