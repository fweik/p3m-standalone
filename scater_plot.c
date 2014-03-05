#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "q.h"
#include "find_error.h"

#define FLOAT_TYPE double;

int main(int argc, char **argv) {
  double err = atof(argv[1]);
  int meth = atoi(argv[2]);

  const int N = 1000;
  const double box = 10.0;
  double factor = 2.0*N*sqrt(1.0/(double)(N)) / ( box*box );
  factor *= factor;
  double alphaL = 5.0;
  

  for( int i = 0; i < n_meshes; i++)
    for( int j = 0; j < n_caos; j++)
      {
	if(factor*p3m_find_error(alphaL, meshes[i], caos[j], meth) <= err*err) {
	  printf("%d %d\n", meshes[i], caos[j]);
	}
      }
}
