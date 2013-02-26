#include <string.h>
#include <stdio.h>
#include <math.h>

#include "types.h"
#include "generate_system.h"
#include "tuning.h"

#include "p3m-ik.h"

#include "common.h"
#include "io.h"

int main(int argc, char **argv) {
  system_t *s;
  parameters_t p;
  int start, stop, step;  
  FLOAT_TYPE prec;
  method_t method = method_p3m_ik;
  FLOAT_TYPE time;

  start = atoi(argv[1]);
  stop = atoi(argv[2]);
  step = atoi(argv[3]);
  prec = atof(argv[4]);

  char *filename;

  for(int i = start; i <= stop; i+= step) {
    s = generate_system( 1, i, 10.0, 1.0);

    printf("Tuning for %d particles.\n", i);

    memset(&p, 0, sizeof(parameters_t));
    p.rcut = 3.0;

    time = Tune( &method, s, &p, prec);
    if( time < 0.0) {
      puts("Tuning failed.");
      continue;
    }
      
    Calculate_reference_forces( s, &p );

    filename = malloc( (int)(6+log(i)/log(10))*sizeof(char));

    sprintf( filename, "%d.dat", i);

    Write_system_cuda( s, &p, filename);

    Free_system(s);
  }
}
