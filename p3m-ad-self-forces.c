#include "types.h"
#include "common.h"
#include "p3m-common.h"

#include <math.h>

#define P3M_SELF_BRILLOUIN 1

/* This is an implementation of equation (9) of
   V. Ballenegger et al., Computer Physics Communications 182(2011)
   The directional indices of the coefficent \beta are called p,q,r
   for clearity in the code. The mx, ... are the components of m' in
   the paper.
*/

FLOAT_TYPE P3M_k_space_calc_self_force( system_t *s, parameters_t *p, data_t *d,
					int m1, int m2, int m3, int dir)
{ 
  int nx, ny, nz;
  int mx, my, mz;
  //  int[3] n;
  int mesh = p->mesh;
  int true_nx,true_ny,true_nz;
  double theSumOverK = 0.0, mesh_i = 1./mesh;
  double U,U_shiftx;
  int P3M_BRILLOUIN_LOCAL = P3M_SELF_BRILLOUIN;

  for(nx=0; nx<mesh; nx++) 
    for(ny=0; ny<mesh; ny++) 
      for(nz=0; nz<mesh; nz++) {
	if ((nx!=0) || (ny!=0) || (nz!=0)) {
	  true_nx = d->nshift[nx];
	  true_ny = d->nshift[ny];
	  true_nz = d->nshift[nz];
	  for (mx=-P3M_BRILLOUIN_LOCAL; mx<=P3M_BRILLOUIN_LOCAL; mx++) {
	    for (my=-P3M_BRILLOUIN_LOCAL; my<=P3M_BRILLOUIN_LOCAL; my++) {
	      for (mz=-P3M_BRILLOUIN_LOCAL; mz<=P3M_BRILLOUIN_LOCAL; mz++) {
		if ((nx==0) || (ny==0) || (nz==0))
		  continue;
		U = pow(sinc(mesh_i * (true_nx + mx*mesh))*sinc(mesh_i * (true_ny  + my*mesh))*sinc(mesh_i * (true_nz + mz*mesh)), p->cao);
		U_shiftx = pow(sinc(mesh_i * (true_nx + (mx+m1)*mesh))*sinc(mesh_i * (true_ny + (my+m2)*mesh))*sinc(mesh_i * (true_nz + (mz+m3)*mesh)), p->cao);

		theSumOverK += (d->G_hat[r_ind( nx, ny, nz)]) * U * U_shiftx;

		//		printf("%d %d %d: U %e U_shiftx %e G_hat() %e\n theSumOverK %e\n", nx, ny, nz, U, U_shiftx, d->G_hat[c_ind( nx, ny, nz)],		       theSumOverK);
		if(isnan(theSumOverK)) {
		  printf("%d %d %d is nan\n", m1, m2, m3);
		  exit(0);
		}
	      }
	    }
	  }
	}
      }
  // Je mets un facteur (-1/2) pour tre en accord avec la formule finale (9) de l'article pour le coefficient b_beta^(m)
  switch(dir) {
    case 0:
      return (PI*m1*mesh_i) / SQR(s->length) * theSumOverK;
      break;
    case 1:
      return (PI*m2*mesh_i) / SQR(s->length) * theSumOverK;
      break;
    case 2:
      return (PI*m3*mesh_i) / SQR(s->length) * theSumOverK;
      break;
  }
  return 0.0;
}

void Init_self_forces( system_t *s, parameters_t *p, data_t *d ) {
  int m1,m2,m3,dim, indp;

  d->self_force_corrections = Init_array(my_power(1+2*P3M_SELF_BRILLOUIN, 3), 3*sizeof(FLOAT_TYPE));

#pragma omp parallel for collapse(4) private(indp)
  for(m1 = -P3M_SELF_BRILLOUIN; m1<=P3M_SELF_BRILLOUIN; m1++)
    for(m2 = -P3M_SELF_BRILLOUIN; m2<=P3M_SELF_BRILLOUIN; m2++)
      for(m3 = -P3M_SELF_BRILLOUIN; m3<=P3M_SELF_BRILLOUIN; m3++) {
	for(dim = 0; dim<3; dim++) {
	  indp = 3*(9*(m1+P3M_SELF_BRILLOUIN)*P3M_SELF_BRILLOUIN*P3M_SELF_BRILLOUIN + 3*(P3M_SELF_BRILLOUIN+m2)*P3M_SELF_BRILLOUIN + (m3+P3M_SELF_BRILLOUIN)) + dim;
	  d->self_force_corrections[indp] = P3M_k_space_calc_self_force( s, p, d, m1, m2, m3, dim);
	/* printf("b(%d, %d, %d) = (%e, %e, %e)\n", i[0], i[1], i[2], d->self_force_corrections[ind-3], */
	/*        d->self_force_corrections[ind-2], d->self_force_corrections[ind-1]); */
	}
      }
#pragma omp barrier
}

void Substract_self_forces( system_t *s, parameters_t *p, data_t *d, forces_t *f ) {
  
  int m[3], dir;
  int id, ind;
  FLOAT_TYPE sin_term;
  FLOAT_TYPE h = s->length / p->mesh;

  for(id=0;id<s->nparticles;id++) {
    ind = 0;
    for(m[0] = -P3M_SELF_BRILLOUIN; m[0]<=P3M_SELF_BRILLOUIN; m[0]++)
      for(m[1] = -P3M_SELF_BRILLOUIN; m[1]<=P3M_SELF_BRILLOUIN; m[1]++)
	for(m[2] = -P3M_SELF_BRILLOUIN; m[2]<=P3M_SELF_BRILLOUIN; m[2]++) {
	  sin_term = SIN(2*PI*(m[0] * s->p->x[id]/h +
			       m[1] * s->p->y[id]/h +
			       m[2] * s->p->z[id]/h)); 
	  for(dir = 0; dir<3; dir++) {
	    f->f->fields[dir][id] -= SQR(s->q[id]) * d->self_force_corrections[ind++] * sin_term;
	  }
	}
  }
}

/* Internal functions */


