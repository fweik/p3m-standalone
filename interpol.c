#include <stdio.h>
#include <stdlib.h>

#include "types.h"
#include "interpol.h"
#include "common.h"

interpolation_t *Init_interpolation(int ip, int derivatives)
{
  FLOAT_TYPE dInterpol=(FLOAT_TYPE)MaxInterpol;
  FLOAT_TYPE x;
  long   i,j;
  
  interpolation_t *ret;
  ret = Init_array( 1, sizeof(interpolation_t));

  ret->interpol = Init_array( ip + 2, sizeof(FLOAT_TYPE *));

  ret->interpol_d = NULL;
  for( i = 0; i <= ip; i++) {
    ret->interpol[i] = Init_array( 2* MaxInterpol + 1, sizeof(FLOAT_TYPE) );
  }
  ret->interpol[ip+1] = NULL;

  if(derivatives) {
    ret->interpol_d = Init_array( ip + 2, sizeof(FLOAT_TYPE *));
    for( i = 0; i <= ip; i++) {
      ret->interpol_d[i] = Init_array( 2* MaxInterpol + 1, sizeof(FLOAT_TYPE) );
    }
    ret->interpol_d[ip+1] = NULL;
  }

  for(j=0;j<=ip;j++)
    for (i=-MaxInterpol; i<=MaxInterpol; i++)
      {
	x=i/(2.0*dInterpol);
	ret->interpol[j][i+MaxInterpol] = caf_bspline(j, x, ip+1);
	if(derivatives)
	  ret->interpol_d[j][i+MaxInterpol] = caf_bspline_d(j, x, ip+1);
      }
  ret->U_hat = caf_bspline_k;
  return ret;
}

//@TODO: remove memleak
void Free_interpolation(interpolation_t *i) {
  fftw_free(i);
}
