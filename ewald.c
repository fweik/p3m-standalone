/* Computes the Ewald sum of a system of point charges. */
/* Time-stamp: <2007-07-31 11:57 olenz> */

 
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
#define kmax 30

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

const method_t method_ewald = { METHOD_EWALD, "Ewald summation.", 
				METHOD_FLAG_G_hat, 
				&Ewald_init, &Ewald_compute_influence_function, &Ewald_k_space, &Ewald_estimate_error };


/*----------------------------------------------------------------------*/
/* GLOBAL VARIABLES */
/*----------------------------------------------------------------------*/

/* Some global variables, internal to ewald.c. */

static int     kmax2; /* reciprocal space cutoff and its square */

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

  fak1 = 2.0/SQR(s->length);
  fak2 = SQR(PI/(p->alpha*s->length));
  
  for (nx=0; nx <= kmax; nx++)
    for (ny=0; ny <= kmax; ny++)
      for (nz=0; nz <= kmax; nz++) {
	n_sqr = SQR(nx) + SQR(ny) + SQR(nz);
	if ((nx==0 && ny==0 && nz==0) || n_sqr > kmax2)
	  d->G_hat[r_ind(nx,ny,nz)] = 0.0;
	else {
	  d->G_hat[r_ind(nx,ny,nz)] = fak1/n_sqr * exp(-fak2*n_sqr);
	}
      }
}  

data_t *Ewald_init(system_t *s, parameters_t *p)
{
  p->mesh = kmax+1;


  data_t *d = Init_data( &method_ewald, s, p );
  kmax2     = SQR(kmax);

  d->mesh = kmax+1;
  
  return d;
}

void Ewald_k_space(system_t *s, parameters_t *p, data_t *d, forces_t *f)
{
  int    i;
  int    nx, ny, nz;
  FLOAT_TYPE kr;
  FLOAT_TYPE rhohat_re, rhohat_im, ghat;
  FLOAT_TYPE force_factor;
  FLOAT_TYPE Leni = 1.0/s->length;
  
  for (nx=-kmax; nx<=kmax; nx++)
    for (ny=-kmax; ny<=kmax; ny++)
      for (nz=-kmax; nz<=kmax; nz++)
	if (nx*nx + ny*ny + nz*nz <= kmax2) {
	  /* compute rhohat */
	  rhohat_re = 0.0;
	  rhohat_im = 0.0;
          ghat = d->G_hat[r_ind(abs(nx), abs(ny), abs(nz))];          
          
#pragma omp parallel for private(kr) reduction( + : rhohat_re) reduction( + : rhohat_im)
	  for (i=0; i<s->nparticles; i++) {
            kr = 2.0*PI*Leni*(nx*(s->p->x[i]) + ny*(s->p->y[i]) + nz*(s->p->z[i]));
	    rhohat_re += s->q[i] * cos(kr);
	    rhohat_im += s->q[i] * -sin(kr);
          }



	  /* compute forces */
#pragma omp parallel for private(kr, force_factor)
	  for (i=0; i<s->nparticles; i++) {
	    kr = 2.0*PI*Leni*(nx*(s->p->x[i]) + ny*(s->p->y[i]) + nz*(s->p->z[i]));
	     
            force_factor = s->q[i] * ghat 
              * (rhohat_re*sin(kr) + rhohat_im*cos(kr));
	      
            f->f_k->x[i] += nx * force_factor;
            f->f_k->y[i] += ny * force_factor;
            f->f_k->z[i] += nz * force_factor;
	    
	  }
	}

  return;
}

FLOAT_TYPE compute_error_estimate_r(system_t *s, parameters_t *p, FLOAT_TYPE alpha) {
  FLOAT_TYPE res;
  FLOAT_TYPE rmax2 = p->rcut*p->rcut;
  /* Kolafa-Perram, eq. 16 */
  res = s->q2*sqrt(p->rcut/(2.0*s->length*s->length*s->length)) * exp(-SQR(alpha)*rmax2) / (SQR(alpha)*rmax2);
  
  return res;
}

FLOAT_TYPE compute_error_estimate_k(system_t *s, parameters_t *p, FLOAT_TYPE alpha) {
  /* compute the k space part of the error estimate */
  FLOAT_TYPE res, Leni = 1.0/s->length;

  /* Kolafa Perram, eq. 31 */
/*   res = Q2 * alpha * my_power(PI, -2.0) * my_power(kmax, -1.5)  */
/*     * exp(-SQR(PI*kmax/(alpha*L))); */

  /* Petersen 1995, eq. 11 */
  res = 2.0 * s->q2 * alpha * Leni * sqrt(1.0/(PI*kmax*s->nparticles))
    * exp(-SQR(PI*kmax/(alpha*s->length)));

  return res;
}

FLOAT_TYPE Ewald_estimate_error(system_t *s, parameters_t *p) {
  return sqrt(SQR(compute_error_estimate_r(s, p, p->alpha)) + SQR(compute_error_estimate_k(s, p, p->alpha)));
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

