#include "statistics.h"
#include "types.h"
#include "common.h"

FLOAT_TYPE *radial_charge_distribution(FLOAT_TYPE r_min, FLOAT_TYPE r_max, int bins, system_t *s) {
  FLOAT_TYPE *rdf = Init_array( bins, 2*sizeof(FLOAT_TYPE));
  FLOAT_TYPE dr = (r_max - r_min) / bins, dx, dy, dz, r2, r, r_min2, r_max2;
  FLOAT_TYPE lengthi = 1.0/s->length;
  int n=0, bin;
  FLOAT_TYPE r_in, r_out;
  FLOAT_TYPE bin_volume;

  memset(rdf, 0, 2*bins*sizeof(FLOAT_TYPE));

  r_min2 = SQR(r_min);
  r_max2 = SQR(r_max);

  for(int i=0;i<s->nparticles;i++) {
    for(int j=0;j<s->nparticles;j++) {
      if(i==j)
	continue;
      if(s->q[i] == s->q[j])
	continue;
      dx = s->p->x[i] - s->p->x[j];
      dx -= ROUND(dx*lengthi)*s->length;
      dy = s->p->x[i] - s->p->y[j];
      dy -= ROUND(dy*lengthi)*s->length;
      dz = s->p->x[i] - s->p->z[j];
      dz -= ROUND(dz*lengthi)*s->length;

      r2 = SQR(dx) + SQR(dy) + SQR(dz);

      if((r2 < r_min2) || (r2 > r_max2))
	continue;

      r = SQRT(r2);

      bin = (int)((r - r_min)/dr);

      rdf[2*bin+1]++;
      n++;
    }

  }

  for(int i=0;i<bins;i++) {
    r_in = i*dr + r_min;
    r_out = r_in + dr;

    bin_volume = 4.0/3.0 * PI * (r_out *r_out*r_out - r_in * r_in * r_in );
    rdf[2*i] = r_in;
    rdf[2*i+1] /= n * bin_volume;
  }
  return rdf;
}
