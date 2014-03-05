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

#include "common.h"
#include "types.h"

#include "timings.h"

#include "realpart.h"

#include "generate_system.h"

#include "ewald.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  int i;
  system_t *s;
  forces_t *f;
  parameters_t p;
  data_t d;
  double error;
  
  p.rcut = 3.0;
  
  double time_nlist, time_n2;
  
  for(i = 1000;i <= 10000;i+=1000) {
    s = generate_system( SYSTEM_RANDOM, i, 10.0*pow( i, 1.0/3.0), 1.0);
    f = Init_forces(i);
    
    d.neighbor_list = malloc(i*sizeof(neighbor_list_t));

    Init_neighborlist( s, &p, &d );
    
    start_timer();
    Realpart_neighborlist( s, &p, &d, f );
    time_nlist = stop_timer();

    error = method_ewald.Error( s, &p );

    //   start_timer();
    //Realteil( s, &p, f );
    //time_n2 = stop_timer();
    
    printf("%d %lf %lf\n", i, time_nlist);
    
    free(d.neighbor_list);
    Free_forces(f);
    Free_system(s);
    
  }
  
}
