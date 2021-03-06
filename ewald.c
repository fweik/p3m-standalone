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

/* Computes the Ewald sum of a system of point charges. */
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "common.h"
#include "p3m-common.h"

#include "ewald.h"

/*----------------------------------------------------------------------*/
/* CONSTANTS */
/*----------------------------------------------------------------------*/
/* Maximal k-vector allowed in Ewald_k_space. 
   Defines the size of the influence function array.
*/

/* Square of the maximal k-vector used for the computation of the
   Kolafa-Perram diagonal term. */
#define KMAX_KPDIAG 50
#define KMAX_KPDIAG_SQR (KMAX_KPDIAG*KMAX_KPDIAG)

/* maximal precision (used in cutoff) */
#define PREC 1.0e-30

/* constant to use when alpha is to be optimized */
#define ALPHA_OPT 0.0

/* precision of alpha */
#define ALPHA_OPT_PREC 1.0e-10

// Method declaration

const method_t method_ewald = { METHOD_EWALD, "Ewald summation.", "ewald",
				METHOD_FLAG_G_hat, 
				&Ewald_init, &Ewald_compute_influence_function, &Ewald_k_space, &Ewald_estimate_error, &Ewald_error_k };

static FLOAT_TYPE compute_error_estimate_k(system_t *s, parameters_t *p, FLOAT_TYPE alpha);

FLOAT_TYPE Ewald_error_k( system_t *s, parameters_t *p ) {
  return compute_error_estimate_k(s, p, p->alpha);
}

FLOAT_TYPE Ewald_self_energy( system_t *s, parameters_t *p ) {
  return - (p->alpha / sqrt(PI)) * s->q2;
}

/*----------------------------------------------------------------------*/
/* GLOBAL VARIABLES */
/*----------------------------------------------------------------------*/

/* Some global variables, internal to ewald.c. */

/*----------------------------------------------------------------------*/
/* HELPER FUNCTIONS */
/*----------------------------------------------------------------------*/
/* Precomputes the influence function 

     2.0/L^2 * exp(-(PI*n/(alpha*L))^2)/n^2 

   as a function of lattice vector n (NOT k=2*PI*n/L). 
   This is stored in the array

     Ghat[Maxkmax+1][Maxkmax+1][Maxkmax+1]

   For symmetry reasons only one octant is actually needed. 
*/
void Ewald_compute_influence_function(system_t *s, parameters_t *p, data_t *d)
{
  
  int    nx,ny,nz;
  FLOAT_TYPE n_sqr,fak1,fak2;
  const int kmax = p->mesh - 1;
  const int kmax2 = kmax*kmax;
  
  fak1 = 2.0/SQR(s->length);
  fak2 = SQR(PI/(p->alpha*s->length));

#ifdef _OPENMP
#pragma omp parallel for collapse(2) private(n_sqr)
#endif 
  for (nx=0; nx <= kmax; nx++)
    for (ny=0; ny <= kmax; ny++)
      for (nz=0; nz <= kmax; nz++) {
	n_sqr = SQR(nx) + SQR(ny) + SQR(nz);
	if ((nx==0 && ny==0 && nz==0) || (n_sqr > kmax2)) {
  //  printf("%d %d %d, kmax %d, kmax2 %d\n", nx, ny, nz, kmax, kmax2);
	  d->G_hat[r_ind(nx,ny,nz)] = 0.0;
        } else {
	  d->G_hat[r_ind(nx,ny,nz)] = fak1/n_sqr * EXP(-fak2*n_sqr);
	}
      }
}  

data_t *Ewald_init(system_t *s, parameters_t *p)
{
  p->tuning = 0;
  p->mesh = p->mesh + 1; 
  data_t *d = Init_data( &method_ewald, s, p );
  d->mesh = p->mesh;
  p->mesh--;

  return d;
}

void Ewald_k_space(system_t *s, parameters_t *p, data_t *d, forces_t *f)
{
  int    i;
  int    nx, ny, nz;
  FLOAT_TYPE kr, energy=0.0;
  FLOAT_TYPE rhohat_re=0, rhohat_im=0, ghat=0;
  FLOAT_TYPE force_factor;
  FLOAT_TYPE Leni = 1.0/s->length;
  const int kmax = p->mesh;
  const int kmax2 = kmax*kmax;
  
  for (nx=-kmax; nx<=kmax; nx++)
    for (ny=-kmax; ny<=kmax; ny++)
      for (nz=-kmax; nz<=kmax; nz++)
	if (nx*nx + ny*ny + nz*nz <= kmax2) {
	  /* compute rhohat */
	  rhohat_re = 0.0;
	  rhohat_im = 0.0;
          ghat = d->G_hat[r_ind(abs(nx), abs(ny), abs(nz))];          
#ifdef _OPENMP          
#pragma omp parallel for private(kr) reduction( + : rhohat_re ) reduction( + : rhohat_im )
#endif
	  for (i=0; i<s->nparticles; i++) {
            kr = 2.0*PI*Leni*(nx*(s->p->x[i]) + ny*(s->p->y[i]) + nz*(s->p->z[i]));
	    rhohat_re += s->q[i] * COS(kr);
	    rhohat_im += s->q[i] * -SIN(kr);
          }
#ifdef _OPENMP
#pragma omp barrier

#pragma omp parallel for private(kr, force_factor)
#endif
	  for (i=0; i<s->nparticles; i++) {
	    kr = 2.0*PI*Leni*(nx*(s->p->x[i]) + ny*(s->p->y[i]) + nz*(s->p->z[i]));
	     
            force_factor = s->q[i] * ghat 
              * (rhohat_re*SIN(kr) + rhohat_im*COS(kr));

            f->f_k->x[i] += nx * force_factor;
            f->f_k->y[i] += ny * force_factor;
            f->f_k->z[i] += nz * force_factor;
	  }
	}
#ifdef _OPENMP
#pragma omp barrier
#endif
  /* compute energy */
  energy += ghat*(SQR(rhohat_re) + SQR(rhohat_im));

  s->energy +=(0.5*s->length /(2.0*PI)) * energy;
  s->energy += Ewald_self_energy(s, p);
}

FLOAT_TYPE compute_error_estimate_r(system_t *s, parameters_t *p, FLOAT_TYPE alpha) {
  FLOAT_TYPE res;
  FLOAT_TYPE rmax2 = p->rcut*p->rcut;
  /* Kolafa-Perram, eq. 16 */
  res = s->q2*SQRT(p->rcut/(2.0*s->length*s->length*s->length)) * exp(-SQR(alpha)*rmax2) / (SQR(alpha)*rmax2);
  
  return res;
}

static FLOAT_TYPE compute_error_estimate_k(system_t *s, parameters_t *p, FLOAT_TYPE alpha) {
  /* compute the k space part of the error estimate */
  FLOAT_TYPE res, Leni = 1.0/s->length;
  const int kmax = p->mesh;

  /* Kolafa Perram, eq. 31 */
/*   res = Q2 * alpha * my_power(PI, -2.0) * my_power(kmax, -1.5)  */
/*     * exp(-SQR(PI*kmax/(alpha*L))); */

  /* Petersen 1995, eq. 11 */
  res = 2.0 * s->q2 * alpha * Leni * SQRT(1.0/(PI*kmax*s->nparticles))
    * EXP(-SQR(PI*kmax/(alpha*s->length)));

  return res;
}

FLOAT_TYPE Ewald_estimate_error(system_t *s, parameters_t *p) {
  return SQRT(SQR(compute_error_estimate_r(s, p, p->alpha)) + SQR(compute_error_estimate_k(s, p, p->alpha)));
}

FLOAT_TYPE Ewald_compute_optimal_alpha(system_t *s, parameters_t *p) {
  /* use bisectional method to get optimal alpha value */
  FLOAT_TYPE alpha_low, f_low;
  FLOAT_TYPE alpha_high, f_high;
  FLOAT_TYPE alpha_guess, f_guess;

  alpha_low = 0.01;
  alpha_high = 10.0;

  f_low = compute_error_estimate_r(s, p, alpha_low) - compute_error_estimate_k(s,p,alpha_low);
  f_high = compute_error_estimate_r(s,p,alpha_high) - compute_error_estimate_k(s,p,alpha_high);

  if (f_low*f_high > 0.0) {
    fprintf(stderr, "Error: Could not init method to find optimal alpha!\n");
    exit(1);
  }

  do {
    alpha_guess = 0.5 *(alpha_low + alpha_high);
    f_guess = compute_error_estimate_r(s, p, alpha_guess) - compute_error_estimate_k(s, p, alpha_guess);
    if (f_low*f_guess < 0.0) alpha_high = alpha_guess;
    else alpha_low = alpha_guess;
  } while (fabs(alpha_low-alpha_high) > ALPHA_OPT_PREC);

  return 0.5 *(alpha_low + alpha_high);
}

