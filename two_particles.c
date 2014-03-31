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
#include "io.h"
#include "p3m-ik.h"
#include "ewald.h"

#include <math.h>

int main(int argc, char *argv[]) {
  FLOAT_TYPE box = 10.0;
  FLOAT_TYPE min_dist = 0.1;
  FLOAT_TYPE max_dist = 8.0;
  FLOAT_TYPE dist_step = 0.01;
  FLOAT_TYPE dist;
  method_t method = method_ewald; 
  parameters_t parameters;
  system_t *system = Init_system(2);
  data_t *d;
  forces_t *f = Init_forces(2);

  parameters.rcut = 1.0;
  parameters.cao = 7;
  parameters.cao3 = parameters.cao*parameters.cao*parameters.cao;
  parameters.ip = parameters.cao - 1;
  parameters.alpha = 0.5;
  parameters.mesh = 32;

  system->q[0] = 1.0;
  system->q[1] = -1.0;
  system->length = 100.0;
  system->q2 = 2.;
  
  system->p->x[0] = system->p->x[1] = 0.5*box; 
  system->p->y[0] = system->p->y[1] = 0.5*box; 
  system->p->z[0] = 0.5*box - 0.5*min_dist;
  system->p->z[1] = 0.5*box + 0.5*min_dist; 
 
  d = method.Init( system, &parameters );
  method.Influence_function( system, &parameters, d );

  while( (dist = FLOAT_ABS(system->p->z[0] - system->p->z[1])) <= max_dist) {
    Calculate_forces( &method, system, &parameters, d, f);
    /* Calculate_reference_forces( system, &parameters ); */

    printf("%lf %e %e\n", dist, f->f->z[0],system->reference->f->z[0]);

    system->p->z[0] -= 0.5*dist_step;
    system->p->z[1] += 0.5*dist_step;    
  }

}
