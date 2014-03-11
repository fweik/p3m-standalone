/**    Copyright (C) 2011,2012,2013 Florian Weik <fweik@icp.uni-stuttgart.de>

       This program is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

       This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this program.  If not, see <http://www.gnu.org/licenses/>. **/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>

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
  FLOAT_TYPE t;
  int m_id = -1;

  start = atoi(argv[1]);
  stop = atoi(argv[2]);
  step = atoi(argv[3]);
  prec = atof(argv[4]);
  rcut = atof(argv[5]);
  density = atof(argv[6]);
  charge = atof(argv[7]);

  if(argc == 9) {
    m_id = atoi(argv[8]);
  }

  char hostname[255];

  gethostname(hostname, sizeof(hostname));

  printf("Running on '%s'\n", hostname);

  for(int i = 0; i < 6; i++) {
    if((m_id != -1) && (i != m_id) )
      continue;
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
    forces_t *forces = Init_forces( s->nparticles );

    for(int j = 0; j < 6; j++) {
      if((m_id != -1) && (j != m_id))
	continue;

      memset(&p, 0, sizeof(parameters_t));
      p.rcut = rcut;
      p.tuning = 1;

      runtime_stat_t timing;

      printf("\t%s:\n", methods[j].method_name);

      t = MPI_Wtime();
      timing = Tune( methods+j, s, &p, prec);
      t = MPI_Wtime() - t;
      if( timing.t.avg < 0.0) {
	printf("\t\tTuning failed.\n");
	continue;
      }
      printf("\t\tmesh %d cao %d time %lf (t_c %e t_f %e t_g %e) (tuning time %lf)\n", p.mesh, p.cao, timing.t.avg, timing.t_c.avg, timing.t_f.avg, timing.t_g.avg, t);

      double tt;
      runtime_t mt;
      if(1) {

	data_t *d = methods[j].Init( s, &p );

	tt = MPI_Wtime();
	methods[j].Kspace_force( s, &p, d, forces );
	tt = MPI_Wtime() - tt;

	printf("\t\tMeasured total time %e (t_c %e t_f %e t_g %e)\n", tt, d->runtime.t_c, d->runtime.t_f, d->runtime.t_g);

	mt = d->runtime;

	mt.t = mt.t_c + mt.t_g + mt.t_f;

	Free_data(d);
      }
    
      fprintf(f[j], "%d %d %d %lf %e %e %e %e %e %e %e %e %e %e %e %e\n", i, p.mesh, p.cao, p.alpha, tt, 
	      timing.t.avg, timing.t.min, timing.t.max,
	      timing.t_c.avg, timing.t_g.avg, timing.t_f.avg,
	      mt.t, mt.t_c, mt.t_g, mt.t_f, p.precision);
      fflush(f[j]);
    }
    
    Free_forces(forces);
    Free_system(s);
  }

  fclose(sys_params);
  for(int i = 0; i < 6; i++) {
    if((m_id != -1) && (i != m_id))
      continue;
    fclose(f[i]);
  }
}
