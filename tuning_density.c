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

#include "p3m-ik-i.h"
#include "p3m-ad-i.h"
#include "p3m-ik-real.h"
#include "p3m-ad-real.h"
#include "ewald.h"
#include "io.h"

#define REFERENCE_SAFTY_FACTOR 1e-3

int main(int argc, char **argv) {
  system_t *s;
  parameters_t p;
  FLOAT_TYPE start=0.1, stop=10.0, step=0.1;  
  FLOAT_TYPE prec=1e-3;
  const method_t methods[] = { method_p3m_ik_r };
  const int n_methods = sizeof(methods)/sizeof(method_t);
  FILE *f[n_methods];
  FLOAT_TYPE box;
  FLOAT_TYPE rcut=3.0;
  FLOAT_TYPE charge = 1.0;
  FLOAT_TYPE t;
  int m_id = -1;

  int npart = 10;

  for(double density = start; density <= stop; density+= step) {
    printf("Tuning for density %e.\n", density);
    box = pow((double)(npart)/density, 0.3333);
    printf("density %lf (%lf), box %lf npart %d \n", FLOAT_CAST (double)(npart)/(box*box*box), FLOAT_CAST density, FLOAT_CAST box, npart);    

    s = generate_system( SYSTEM_RANDOM, npart, box, charge);
    forces_t *forces = Init_forces( s->nparticles );
    forces_t *forces_ewald = Init_forces( s->nparticles );
    data_t *d_ewald = NULL;

    for (int i=0; i<3; i++ ) {
        memset ( forces->f->fields[i]  , 0, s->nparticles*sizeof ( FLOAT_TYPE ) );
        memset ( forces->f_k->fields[i], 0, s->nparticles*sizeof ( FLOAT_TYPE ) );
        memset ( forces->f_r->fields[i], 0, s->nparticles*sizeof ( FLOAT_TYPE ) );
        memset ( forces_ewald->f->fields[i]  , 0, s->nparticles*sizeof ( FLOAT_TYPE ) );
        memset ( forces_ewald->f_k->fields[i], 0, s->nparticles*sizeof ( FLOAT_TYPE ) );
        memset ( forces_ewald->f_r->fields[i], 0, s->nparticles*sizeof ( FLOAT_TYPE ) );
    }


    if(1) {
    FLOAT_TYPE V = pow( box, 3);
    double alpha = SQRT(-LOG((prec*SQRT(s->nparticles*rcut*V))/(2*SQRT(2)*s->q2)))/rcut;
    parameters_t ewald_parameters;
    ewald_parameters.rcut = rcut;
    ewald_parameters.alpha = alpha;
    ewald_parameters.tuning = 0;
    ewald_parameters.mesh = 4;

    FLOAT_TYPE ref_err;
    while((ref_err = Ewald_error_k( s, &ewald_parameters )) > REFERENCE_SAFTY_FACTOR*prec) {
      ewald_parameters.mesh += 2;
    }

    d_ewald = method_ewald.Init( s, &ewald_parameters);

    method_ewald.Kspace_force( s, &ewald_parameters, d_ewald, forces_ewald);
    
    Free_data(d_ewald);
    }

    for(int j = 0; j < n_methods; j++) {
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
      printf("\t\tmesh %d cao %d alpha %e time %lf prec %e (t_c %e t_f %e t_g %e) (tuning time %lf)\n", p.mesh, p.cao, p.alpha, timing.t.avg, p.precision, timing.t_c.avg, timing.t_f.avg, timing.t_g.avg, t);

      double tt;
      runtime_t mt;
      if(1) {
	p.tuning = 0;
	data_t *d = methods[j].Init( s, &p );

	for(int i = 0; i < 3; i++)
	  memset(forces->f_k->fields[i], 0, s->nparticles*sizeof(FLOAT_TYPE));

	tt = MPI_Wtime();
	methods[j].Kspace_force( s, &p, d, forces );
	tt = MPI_Wtime() - tt;

	printf("\t\tMeasured total time %e (t_c %e t_f %e t_g %e)\n", tt, d->runtime.t_c, d->runtime.t_f, d->runtime.t_g);

	mt = d->runtime;

	mt.t = mt.t_c + mt.t_g + mt.t_f;

	double actual_error = 0.0;
 
	for(int i = 0; i < s->nparticles; i++)
	  for(int j = 0; j < 3; j++) {
	    actual_error += SQR(forces->f_k->fields[j][i] - forces_ewald->f_k->fields[j][i]);  
	    // printf("id %d dim %d method_force %e ewald_force %e \n", i, j, forces->f_k->fields[j][i], forces_ewald->f_k->fields[j][i]);
	  }

	printf("actual_error %e\n", SQRT(actual_error / s->nparticles));
	Free_data(d);
      }
    
      printf("out %d %e %d %d %lf %e %e %e %e %e %e %e %e %e %e %e %e\n", npart, density, p.mesh, p.cao, p.alpha, tt, 
	      timing.t.avg, timing.t.min, timing.t.max,
	      timing.t_c.avg, timing.t_g.avg, timing.t_f.avg,
	      mt.t, mt.t_c, mt.t_g, mt.t_f, p.precision);
    }
    
    Free_forces(forces);
    Free_forces(forces_ewald);
    Free_system(s);
  }

  for(int i = 0; i < n_methods; i++) {
    if((m_id != -1) && (i != m_id))
      continue;
    fclose(f[i]);
  }
}
