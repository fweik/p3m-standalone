#include "p3m-ik-cuda.h"
#include "cuda/p3m-ik-cuda.h"
#include "p3m-ik.h"

#include <stdlib.h>

const method_t method_p3m_ik_cuda = { METHOD_P3M_ik_cuda, "P3M with ik differentiation, not intelaced, gpu.",
				      0, &Init_ik_cuda, NULL, &P3M_ik_cuda, &Error_ik, &Error_ik_k };


data_t *Init_ik_cuda( system_t *s, parameters_t *p ) {
  data_t *d = malloc(sizeof(data_t));
  p3m_cuda_data_t *cd = malloc(sizeof(p3m_cuda_data_t));
  double *pos = malloc(3*s->nparticles*sizeof(double));

  d->G_hat = NULL;
  d->Qmesh = NULL;
  d->Fmesh = NULL;
  d->nshift = NULL;
  d->Dn = NULL;
  for(int i = 0; i<2; i++) {
    d->dQdx[i] = NULL;
    d->dQdy[i] = NULL;
    d->dQdz[i] = NULL;
  }

  d->inter = NULL;
  d->cf[0] = d->cf[1] = NULL;
  d->ca_ind[0] = d->ca_ind[1] = NULL;
  d->forward_plans = 0;
  d->backward_plans = 0;

  for(int i = 0; i < s->nparticles; i++) {
    pos[3*i + 0] = s->p->x[i];
    pos[3*i + 1] = s->p->y[i];
    pos[3*i + 2] = s->p->z[i];    
  }

  d->method_data = cd;

  cd->n = s->nparticles;
  cd->box = s->length;
  cd->pos = pos;
  cd->q = s->q;

  cd->alpha = p->alpha;
  cd->cao = p->cao;
  cd->mesh = p->mesh;

  cd->f = NULL;

  printf("Calling p3m_ik_cuda_init for %d particles, cao %d, mesh %d\n", cd->n, cd->cao, cd->mesh);

  p3m_ik_cuda_init( cd );

  return d;
}

void P3M_ik_cuda( system_t *s, parameters_t *p, data_t *d, forces_t *f) {
  p3m_cuda_data_t *cd = (p3m_cuda_data_t *) d->method_data;

  cd->f = f->f_k->fields;

  p3m_ik_cuda( cd );
}


