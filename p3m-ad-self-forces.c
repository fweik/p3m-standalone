/**    Copyright (C) 2011,2012,2013 Florian Weik <fweik@icp.uni-stuttgart.de>

       This program is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

       This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this program.  If not, see <http://www.gnu.org/licenses/>. **/

#include "types.h"
#include "common.h"
#include "p3m-common.h"
#include "p3m-ad-self-forces.h"

#include <math.h>

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
  FLOAT_TYPE theSumOverK = 0.0, mesh_i = 1./mesh;
  FLOAT_TYPE U,U_shiftx;
  int P3M_BRILLOUIN_LOCAL = P3M_SELF_BRILLOUIN;
  FLOAT_TYPE G_hat;

  SF_TRACE(puts("P3M_k_space_calc_self_force()"));

  for(nx=0; nx<mesh; nx++) 
    for(ny=0; ny<mesh; ny++) 
      for(nz=0; nz<mesh; nz++) {
	G_hat = d->G_hat[r_ind(nx, ny, nz)];
	if(G_hat == 0.0)
	  continue;
	true_nx = d->nshift[nx];
	true_ny = d->nshift[ny];
	true_nz = d->nshift[nz];
	//	printf("n[] = (%d %d %d), true_n[] = (%d %d %d)\n", nx, ny, nz, true_nx, true_ny, true_nz);
	#ifdef _OPENMP
#pragma omp parallel for private(U, U_shiftx) reduction( + : theSumOverK ) collapse(3)
#endif
	for (mx=-P3M_BRILLOUIN_LOCAL; mx<=P3M_BRILLOUIN_LOCAL; mx++) {
	  for (my=-P3M_BRILLOUIN_LOCAL; my<=P3M_BRILLOUIN_LOCAL; my++) {
	    for (mz=-P3M_BRILLOUIN_LOCAL; mz<=P3M_BRILLOUIN_LOCAL; mz++) {
	      U = pow(sinc(mesh_i * (true_nx + mx*mesh))*sinc(mesh_i * (true_ny  + my*mesh))*sinc(mesh_i * (true_nz + mz*mesh)), p->cao);
	      U_shiftx = pow(sinc(mesh_i * (true_nx + (mx+m1)*mesh))*sinc(mesh_i * (true_ny + (my+m2)*mesh))*sinc(mesh_i * (true_nz + (mz+m3)*mesh)), p->cao);

	      theSumOverK += G_hat * U * U_shiftx;

	      /* printf("%d %d %d: U %e U_shiftx %e G_hat() %e\n theSumOverK %e\n", nx, ny, nz, FLOAT_CAST U, U_shiftx, FLOAT_CAST d->G_hat[c_ind( nx, ny, nz)],  FLOAT_CAST theSumOverK); */
	      if(isnan(theSumOverK)) {
		printf("%d %d %d is nan\n", m1, m2, m3);
		exit(0);
	      }
	    }
	  }
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
      }
 
  switch(dir) {
  case 0:
    return (PI*m1*p->mesh) / SQR(SQR(s->length)) * theSumOverK;
    break;
  case 1:
    return (PI*m2*p->mesh) / SQR(SQR(s->length)) * theSumOverK;
    break;
  case 2:
    return (PI*m3*p->mesh) / SQR(SQR(s->length)) * theSumOverK;
    break;
  }
  return 0.0;
}

void Init_self_forces( system_t *s, parameters_t *p, data_t *d ) {
  int m[3],dim, indp=0;

  SF_TRACE(puts("Init_self_forces()"));

  d->self_force_corrections = (FLOAT_TYPE *)Init_array(my_power(1+2*P3M_SELF_BRILLOUIN, 3), 3*sizeof(FLOAT_TYPE));

  for(m[0] = -P3M_SELF_BRILLOUIN; m[0]<=P3M_SELF_BRILLOUIN; m[0]++)
    for(m[1] = -P3M_SELF_BRILLOUIN; m[1]<=P3M_SELF_BRILLOUIN; m[1]++)
      for(m[2] = -P3M_SELF_BRILLOUIN; m[2]<=P3M_SELF_BRILLOUIN; m[2]++) {
	for(dim = 0; dim<3; dim++) {
	  d->self_force_corrections[indp++] = (m[dim] == 0) ? 0.0 : P3M_k_space_calc_self_force( s, p, d, m[0], m[1], m[2], dim);
	}
	/* printf("b(%d, %d, %d) = (%e, %e, %e)\n", m[0], m[1], m[2], FLOAT_CAST d->self_force_corrections[indp-3],  */
	/*        FLOAT_CAST d->self_force_corrections[indp-2], FLOAT_CAST d->self_force_corrections[indp-1]);  */
      }
#ifdef _OPENMP
#pragma omp barrier
#endif
}

void Substract_self_forces( system_t *s, parameters_t *p, data_t *d, forces_t *f ) {
  
  int m[3], dir;
  int id, ind;
  FLOAT_TYPE sin_term;
  FLOAT_TYPE h = s->length / p->mesh;
  //  FLOAT_TYPE f_self[3] = { 0.0, 0.0, 0.0};

  for(id=0;id<s->nparticles;id++) {
    ind = 0;
    for(m[0] = -P3M_SELF_BRILLOUIN; m[0]<=P3M_SELF_BRILLOUIN; m[0]++)
      for(m[1] = -P3M_SELF_BRILLOUIN; m[1]<=P3M_SELF_BRILLOUIN; m[1]++)
	for(m[2] = -P3M_SELF_BRILLOUIN; m[2]<=P3M_SELF_BRILLOUIN; m[2]++) {
	  sin_term = SQR(s->q[id]) * SIN(2*PI*(m[0] * s->p->x[id]/h +
			       m[1] * s->p->y[id]/h +
			       m[2] * s->p->z[id]/h)); 
	  /* printf("Selfforce: m (%d %d %d), sin_term %e, self_force_corrections (%e %e %e)\n", m[0], m[1], m[2], FLOAT_CAST sin_term, FLOAT_CAST d->self_force_corrections[ind], */
	  /* 	 FLOAT_CAST d->self_force_corrections[ind+1], FLOAT_CAST d->self_force_corrections[ind+2]); */

	  for(dir = 0; dir<3; dir++) {
	    f->f->fields[dir][id] -=  d->self_force_corrections[ind++] * sin_term;
	  }
	}
    /* printf("Selfforce: particle %d, force (%e %e %e)\n", id, FLOAT_CAST f_self[0], FLOAT_CAST f_self[1], FLOAT_CAST f_self[2]); */
  }
}

/* Internal functions */


