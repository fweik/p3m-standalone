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

#include "wtime.h"
#include "types.h"
#include "generate_system.h"
#include "charge-assign.h"
#include "p3m-common.h"
#include "p3m-ik.h"
#include "p3m-ad-i.h"

int main(int argc, char **argv) {
  system_t *s;
  parameters_t p;
  data_t *d;
  int start, stop, step;  
  int mesh;
  FLOAT_TYPE box;
  FLOAT_TYPE density;
  FLOAT_TYPE t, total;
  FILE *f;
  double sum = 0.0;

  start = atoi(argv[1]);
  stop = atoi(argv[2]);
  step = atoi(argv[3]);
  density = atof(argv[4]);
  mesh = atoi(argv[5]);

  f = fopen("assignment.dat", "w");

  fprintf(f, "#particles ");
  
  for(int j = 1; j <= 7; j++) {
    fprintf(f, "assign_charge-%d assign_forces-%d assign_charge_real-%d assign_forces_real-%d assign_charge_real_nostor-%d assign_forces_nostor-%d", j,j,j,j,j,j);
  }
  fprintf(f, "\n");

  for(int i = start; i <= stop; i+= step) {
    box = pow((double)(i)/density, 0.3333);
    s = generate_system( SYSTEM_RANDOM, i, box, 1.0);
    forces_t *forces = Init_forces(s->nparticles);
    
    printf("Tuning for %d particles box %e, dens %e.\n", s->nparticles, box, (double)(i)/(box*box*box));

    fprintf(f, "%d ", i);

    total = wtime();

    for(int j = 1; j <= 7; j++) {
      memset(&p, 0, sizeof(parameters_t));
      p.tuning = 1;
      p.cao = j;
      p.cao3 = j*j*j;
      p.ip = j-1;
      p.mesh = mesh;
      d = Init_data(&method_p3m_ik, s, &p);

      /* t = wtime(); */
      /* assign_charge(s, &p, d, 0); */
      /* t = wtime() - t; */

      /* fprintf(f, "%e ", t); */

      /* t = wtime(); */
      /* assign_forces(1.0, s, &p, d, s->reference, 0); */
      /* t  = wtime() - t; */
      /* fprintf(f, "%e ", t); */

      /* t = wtime(); */
      /* assign_charge_real(s, &p, d); */
      /* t = wtime() - t; */

      t = wtime();
      assign_charge_real_res(s, &p, d);
      t = wtime() - t;

      fprintf(f, "%e ", t);

      /* t = wtime(); */
      /* assign_charge_and_derivatives_real( s, &p, d); */
      /* t = wtime() - t; */

      /* fprintf(f, "%e ", t); */

      /* t = wtime(); */
      /* assign_charge_and_derivatives_real_res( s, &p, d); */
      /* t = wtime() - t; */

      /* fprintf(f, "%e ", t); */

      t = wtime();
      assign_charge_real_nostor(s, &p, d);
      t = wtime() - t;

      fprintf(f, "%e ", t);

      /* t = wtime(); */
      /* assign_charge_real_nostor_res(s, &p, d); */
      /* t = wtime() - t; */

      /* fprintf(f, "%e ", t); */

      /* t = wtime(); */
      /* assign_charge_real_nostor_res_5(s, &p, d); */
      /* t = wtime() - t; */

      /* fprintf(f, "%e ", t); */

      for(int i = 0; i < 2*mesh*mesh*mesh; i++) {
      	d->Fmesh->x[i] = 1.1*d->Qmesh[i];
      	d->Fmesh->y[i] = 1.2*d->Qmesh[i];
      	d->Fmesh->z[i] = 1.3*d->Qmesh[i];
      }
      

      t = wtime();
      assign_forces_real(1.0, s, &p, d, forces);
      t  = wtime() - t;
      fprintf(f, "%e ", t);

      /* t = wtime(); */
      /* assign_charge_real_nostor(s, &p, d); */
      /* t = wtime() - t; */

      /* fprintf(f, "%e ", t); */

      t = wtime();
      assign_forces_real_nostor(1.0, s, &p, d, forces);
      t  = wtime() - t;
      fprintf(f, "%e ", t);

      /* t = wtime(); */
      /* assign_forces_ad(1.0, s, &p, d, s->reference, 0); */
      /* assign_forces_ad(1.0, s, &p, d, s->reference, 1); */
      /* t = wtime() - t; */

      /* fprintf(f, "%e ", t); */

      /* assign_charge_and_derivatives(s, &p, d, 0); */
      /* assign_charge_and_derivatives(s, &p, d, 1); */

      /* t = wtime(); */
      /* assign_forces_interlacing_ad(1.0, s, &p, d, s->reference); */
      /* t  = wtime() - t; */
      /* fprintf(f, "%e ", t); */
  
      Free_data(d);
      fflush(f);
    }

    fprintf(f, "\n");

    total = wtime() - total;

    printf("total %e\n", total);

    /* for(int l = 0; l < 3; l++) */
    /*   for(int k = 0; k < i; k++) */
    /* 	sum += s->reference->f_k->fields[l][k]; */
    Free_forces(forces);
    Free_system(s);      
  }

  printf("sum %e\n", sum);
  fclose(f);
}
