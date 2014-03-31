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

#include "types.h"
#include "generate_system.h"
#include "charge-assign.h"
#include "p3m-common.h"
#include "p3m-ik.h"

int main(int argc, char **argv) {
  system_t *s;
  parameters_t p;
  data_t *d;
  int start, stop, step;  
  FLOAT_TYPE box;
  FLOAT_TYPE density;
  FLOAT_TYPE charge;
  FLOAT_TYPE t;
  char filename[256];
  FILE *f;

  start = atoi(argv[1]);
  stop = atoi(argv[2]);
  step = atoi(argv[3]);
  density = atof(argv[4]);
  charge = atof(argv[5]);

  for(int i = start; i <= stop; i*= step) {
    printf("Tuning for %d particles.\n", i);
    box = pow((double)(i)/density, 0.3333);

    s = generate_system( SYSTEM_RANDOM, i, box, charge);

    sprintf(filename, "assign-%d.dat", i);

    f = fopen(filename, "w");

    for(int j = 1; j < 7; j++) {
      memset(&p, 0, sizeof(parameters_t));
      p.tuning = 1;
      p.cao = j;
      p.cao3 = j*j*j;
      p.ip = j-1;
      p.mesh = 8;
      d = Init_data(&method_p3m_ik, s, &p);
      
      fprintf(f, "%d ", j);

      t = MPI_Wtime();
      assign_charge(s, &p, d, 0);
      t = MPI_Wtime() - t;

      fprintf(f, "%e ", t);

      t = MPI_Wtime();
      assign_forces(1.0, s, &p, d, s->reference, 0);
      t  = MPI_Wtime() - t;
      fprintf(f, "%e ", t);

      t = MPI_Wtime();
      assign_charge_real(s, &p, d);
      t = MPI_Wtime() - t;

      fprintf(f, "%e ", t);

      t = MPI_Wtime();
      assign_forces_real(1.0, s, &p, d, s->reference);
      t  = MPI_Wtime() - t;
      fprintf(f, "%e ", t);
  
      fprintf(f, "\n");

      Free_data(d);
      fflush(f);
    }

    fclose(f);
    Free_system(s);
  }
}
