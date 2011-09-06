#include "common.h"
#include "types.h"

#include "timings.h"

#include "realpart.h"

#include "generate_system.h"

#include <math.h>
#include <stdio.h>

int main(void) {
  int i;
  system_t *s;
  forces_t *f;
  parameters_t p;
  
  p.rcut = 1.0;
  
  double time_nlist, time_n2;
  
  for(i = 10;i <= 1000000;i*=10) {
    s = generate_system( FORM_FACTOR_RANDOM, i, 10.0*pow( i, 1.0/3.0), 1.0);
    f = Init_forces(i);
    
    Init_neighborlist( s, &p );
    
    start_timer();
    Realpart_neighborlist( s, &p, f );
    time_nlist = stop_timer();
    
    start_timer();
    Realteil( s, &p, f );
    time_n2 = stop_timer();
    
    printf("%d %lf %lf\n", i, time_nlist, time_n2);
    
    Free_forces(f);
    Free_system(s);
    
  }
  
}
