#include "dipol.h"

void Dipol(system_t *s)
{
  FLOAT_TYPE SummeX,SummeY,SummeZ;
  FLOAT_TYPE VorFak;

  int i;

  SummeX=SummeY=SummeZ=0.0;
  VorFak=4.0*PI/(3.0*s->length*s->length*s->length);

#pragma omp parallel for reduction ( + : SummeX ) reduction ( + : SummeY ) reduction ( + : SummeZ )
  for (i=0; i<s->nparticles; ++i)
    if (s->q[i]!=0)
      {
	SummeX += s->q[i]*s->p.x[i];
	SummeY += s->q[i]*s->p.y[i];
	SummeZ += s->q[i]*s->p.z[i];
      }

 //   Coulomb-Energie: 

  SummeX *= VorFak;
  SummeY *= VorFak;  
  SummeZ *= VorFak;

 // Kraefte: 
#pragma omp parallel for
  for (i=0; i<s->nparticles; ++i)
    if (s->q[i]!=0)
      {
	s->f.x[i] -= s->q[i]*SummeX;
	s->f.y[i] -= s->q[i]*SummeY;
	s->f.z[i] -= s->q[i]*SummeZ;
      }
}