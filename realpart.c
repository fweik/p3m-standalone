#include <math.h>

#include "realpart.h"

void Realteil(system_t *s, p3m_parameters_t *p)
{
  /* Zwei Teilchennummern: */
  int t1,t2;
  /* Minimum-Image-Abstand: */
  FLOAT_TYPE dx,dy,dz,r;
  /* Staerke der elegktrostatischen Kraefte */
  FLOAT_TYPE fak;
  /* Zur Approximation der Fehlerfunktion */
  FLOAT_TYPE ar,erfc_teil;
  FLOAT_TYPE lengthi = 1.0/s->length;
  const FLOAT_TYPE wupi = 1.77245385090551602729816748334;

#pragma omp parallel for private(dx,dy,dz, r, ar, erfc_teil, fak, t2)
  for (t1=0; t1<s->nparticles-1; t1++)   /* Quick and Dirty N^2 */
    for (t2=t1+1; t2<s->nparticles; t2++)
      {
	dx = s->p.x[t1] - s->p.x[t2]; dx -= round(dx*lengthi)*s->length;
	dy = s->p.y[t1] - s->p.y[t2]; dy -= round(dy*lengthi)*s->length; 
	dz = s->p.z[t1] - s->p.z[t2]; dz -= round(dz*lengthi)*s->length;

	r = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
	if (r<=p->rcut)
	  {
	    ar= p->alpha*r;
	    //t = 1.0 / (1.0 + p*ar);
	    //erfc_teil = t*(a1+t*(a2+t*(a3+t*(a4+t*a5))));
	    erfc_teil = erfc(ar);
	    fak = s->q[t1]*s->q[t2]*
	      (erfc_teil/r+(2*p->alpha/wupi)*exp(-ar*ar))/SQR(r);
	    
	    s->f_r.x[t1] += fak*dx;
	    s->f_r.y[t1] += fak*dy;
	    s->f_r.z[t1] += fak*dz;
	    s->f_r.x[t2] -= fak*dx;
	    s->f_r.y[t2] -= fak*dy;
	    s->f_r.z[t2] -= fak*dz;
	  }	 
      }
}