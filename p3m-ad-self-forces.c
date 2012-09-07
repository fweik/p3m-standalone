#include "types.h"
#include "common.h"
#include "p3m-common.h"

#define P3M_SELF_BRILLOUIN 1

FLOAT_TYPE Aliasing_sums_ad(int NX, int NY, int NZ, system_t *s, parameters_t *p, data_t *);

void Init_self_forces( system_t *s, parameters_t *p, data_t *d ) {

}

void Substract_self_forces( system_t *s, parameters_t *p, data_t *d ) {

}

/* Internal functions */

/* This is an implementation of equation (9) of
   V. Ballenegger et al., Computer Physics Communications 182(2011)
   The directional indices of the coefficent \beta are called p,q,r
   for clearity in the code. The mx, ... are the components of m' in
   the paper.
*/

FLOAT_TYPE Aliasing_sums_ad_self(int nx, int ny, int nz, 
				 int p, int q, int r,
				 system_t *s, parameters_t *par, data_t *d)
{
  int    mx,my,mz;
  FLOAT_TYPE NMX,NMY,NMZ;
  FLOAT_TYPE S1, S2, S3;
  int Mesh = p->mesh;
  FLOAT_TYPE Len = s->length;
  FLOAT_TYPE Leni = 1.0/Len;
  FLOAT_TYPE Res = 0.0;
  FLOAT_TYPE fak1 = 1.0;

  for (MX = -P3M_SELF_BRILLOUIN; MX <= P3M_SELF_BRILLOUIN; MX++) {
    NMX = NX + Mesh*MX;
    S1   = my_power(sinc(fak1*NMX), 2*p->cao); 
    for (MY = -P3M_SELF_BRILLOUIN; MY <= P3M_SELF_BRILLOUIN; MY++) {
      NMY = NY + Mesh*MY;
      S2   = S1*my_power(sinc(fak1*NMY), 2*p->cao);
      for (MZ = -P3M_SELF_BRILLOUIN; MZ <= P3M_SELF_BRILLOUIN; MZ++) {
	NMZ = NZ + Mesh*MZ;
	S3   = S2*my_power(sinc(fak1*NMZ), 2*p->cao);
	
	Res += S3;
      }
    }
  }
  return Res;
}
