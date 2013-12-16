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

#include "types.h"
#include "common.h"
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <string.h>

#define FOREACH(A,B) for(int A = 0; A < B; A++)
#define IND(x,y,z,N) N*N*x+N*y+z

FLOAT_TYPE inhomo_error(system_t *s, parameters_t *p, int bins) {
  int n[3];
  FLOAT_TYPE h = s->length / bins;
  FLOAT_TYPE *rho_mesh = Init_array( 2*bins*bins*bins, sizeof(FLOAT_TYPE));

  FFTW_PLAN forward_plan;
  FFTW_PLAN backward_plan;

  FILE *debug_output;

  memset(rho_mesh, 0, 2*bins*bins*bins*sizeof(FLOAT_TYPE));

  for(int id=0;id<s->nparticles;id++) {
    for(int j=0;j<3;j++) {
      n[j] = (int)FLOOR(s->p->fields[j][id] / h);
    }
    rho_mesh[2*IND(n[0],n[1],n[2],bins)] += s->q[id]*s->q[id];
  }
  
  /* debug_output = fopen( "rho_mesh.dat", "w" ); */

  /* FOREACH(i, bins) { */
  /*   FOREACH(j, bins) { */
  /*     FOREACH(k, bins) { */
  /* 	fprintf(debug_output, "%d %d %d %lf\n", i, j, k, rho_mesh[bins*bins*i + bins*j + k]); */
  /*     } */
  /*   } */
  /* } */
  /* fclose(debug_output); */

  forward_plan = FFTW_PLAN_DFT_3D(bins, bins, bins, (FFTW_COMPLEX *)rho_mesh,(FFTW_COMPLEX *)rho_mesh, FFTW_FORWARD, FFTW_PATIENT);

  FFTW_EXECUTE(forward_plan);

  return 0;
}


