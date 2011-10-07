#include "common.h"
#include "types.h"

#include "timings.h"

#include "realpart.h"

#include "generate_system.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  int i;
  system_t *s;
  forces_t *f;
  parameters_t p;
  data_t d;
  
  p.rcut = 1.0;
  
  double time_nlist, time_n2;
  
  for(i = 1000;i <= 30000;i+=1000) {
    s = generate_system( FORM_FACTOR_RANDOM, i, 10.0*my_power( i, 1.0/3.0), 1.0);
    f = Init_forces(i);
    
    d.neighbor_list = malloc(i*sizeof(neighbor_list_t));

    Init_neighborlist( s, &p, &d );
    
    start_timer();
    Realpart_neighborlist( s, &p, &d, f );
    time_nlist = stop_timer();
    
    start_timer();
    Realteil( s, &p, f );
    time_n2 = stop_timer();
    
    printf("%d %lf %lf\n", i, time_nlist, time_n2);
    
    free(d.neighbor_list);
    Free_forces(f);
    Free_system(s);
    
  }
  
}
