#include <fftw3.h>

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
      dx = s->p->x[i] - s->p->x[j];
      dx -= ROUND(dx*lengthi)*s->length;
      dy = s->p->y[i] - s->p->y[j];
      dy -= ROUND(dy*lengthi)*s->length;
      dz = s->p->z[i] - s->p->z[j];
      dz -= ROUND(dz*lengthi)*s->length;

      r2 = SQR(dx) + SQR(dy) + SQR(dz);

      if((r2 < r_min2) || (r2 > r_max2))
	continue;

      r = SQRT(r2);

      bin = (int)((r - r_min)/dr);

      if((r == 0) || ( bin < 0) || ( bin >= bins)) {
	printf("dist(%d, %d) = %lf, bin %d\n", i, j, r, bin);
	printf("r1 (%lf %lf %lf), r2 (%lf %lf %lf)\n", s->p->x[i], s->p->y[i],s->p->z[i],s->p->x[j], s->p->y[j],s->p->z[j]);
	printf("dx %lf dy %lf dz %lf\n", dz, dy, dz);
	exit(0);
      }

      rdf[2*bin+1]+= s->q[i] * s->q[j];
      n++;
    }

  }
  printf("#rdf[0] %lf\n", rdf[0]);
  printf("#rdf[1] %lf\n", rdf[1]);
  printf("#rdf[2] %lf\n", rdf[2]);
  for(int i=0;i<bins;i++) {
    r_in = i*dr + r_min;
    r_out = r_in + dr;

    bin_volume = 4.0/3.0 * PI * (r_out *r_out*r_out - r_in * r_in * r_in );
    rdf[2*i] = 0.5 * (r_in + r_out);
    rdf[2*i+1] /= n * bin_volume / (s->length *SQR(s->length));
  }
  return rdf;
}

FLOAT_TYPE *rdf_fft( int N, FLOAT_TYPE *rdf) {
  fftw_plan p;
  FLOAT_TYPE *data = Init_array( N, 2*sizeof(FLOAT_TYPE));

  memset(data, 0, 2*N*sizeof(FLOAT_TYPE));

  for(int i=0;i<N;i++)
    data[2*i] = rdf[2*i+1];

  p = fftw_plan_dft_1d( N, (fftw_complex *)data,  (fftw_complex *)data, FFTW_FORWARD, FFTW_PATIENT); 

  FFTW_EXECUTE(p);

  return data;
}
