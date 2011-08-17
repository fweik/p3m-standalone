/* Computes the Ewald sum of a system of point charges. */
/* Time-stamp: <2007-07-31 11:57 olenz> */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

/*----------------------------------------------------------------------*/
/* GLOBAL VARIABLES */
/*----------------------------------------------------------------------*/

/* Some global variables, internal to ewald.c. */

static int     kmax2; /* reciprocal space cutoff and its square */
static FLOAT_TYPE Ghat[kmax+1][kmax+1][kmax+1];

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
void Ewald_compute_influence_function(FLOAT_TYPE alpha)
{
  
  int    nx,ny,nz;
  FLOAT_TYPE n_sqr,fak1,fak2;

  fak1 = 2.0/SQR(Len);
  fak2 = SQR(PI/(alpha*Len));
  
  for (nx=0; nx <= kmax; nx++)
    for (ny=0; ny <= kmax; ny++)
      for (nz=0; nz <= kmax; nz++) {
	n_sqr = SQR(nx) + SQR(ny) + SQR(nz);
	if ((nx==0 && ny==0 && nz==0) || n_sqr > kmax2)
	  Ghat[nx][ny][nz] = 0.0;
	else {
	  Ghat[nx][ny][nz] = fak1/n_sqr * exp(-fak2*n_sqr);
	}
      }
}  

void Ewald_init(int particlenumber)
{
  kmax2     = SQR(kmax);
}

void Ewald_k_space(system_t *s, p3m_parameters_t *p)
{
  int    i;
  int    nx, ny, nz;
  FLOAT_TYPE kr;
  FLOAT_TYPE rhohat_re, rhohat_im, ghat;
  FLOAT_TYPE force_factor;

  for (nx=-kmax; nx<=kmax; nx++)
    for (ny=-kmax; ny<=kmax; ny++)
      for (nz=-kmax; nz<=kmax; nz++)
	if (nx*nx + ny*ny + nz*nz <= kmax2) {
	  /* compute rhohat */
	  rhohat_re = 0.0;
	  rhohat_im = 0.0;
#pragma omp parallel for private(kr) reduction( + : rhohat_re) reduction( + : rhohat_im)
	  for (i=0; i<NP; i++) {
            kr = 2.0*PI*Leni*(nx*s->p.x[i] + ny*s->p.y[i] + nz*s->p.z[i];
	    rhohat_re += s->q[i] * cos(kr);
	    rhohat_im += s->q[i] * -sin(kr);
          }

	  ghat = Ghat[abs(nx)][abs(ny)][abs(nz)];

	  /* compute forces if requested */
#pragma omp parallel for private(kr, force_factor)
	  for (i=0; i<NP; i++) {
	    kr = 2.0*PI*Leni*(nx*xS[i] + ny*yS[i] + nz*zS[i]);
	     
            force_factor = Q[i] * ghat 
              * (rhohat_re*sin(kr) + rhohat_im*cos(kr));
	      
            Fx_K[i] += nx * force_factor;
            Fy_K[i] += ny * force_factor;
            Fz_K[i] += nz * force_factor;
	    
	  }
	}

  return;
}

FLOAT_TYPE compute_error_estimate_r(FLOAT_TYPE alpha, FLOAT_TYPE rmax) {
  FLOAT_TYPE res;
  FLOAT_TYPE rmax2 = rmax*rmax;
  /* Kolafa-Perram, eq. 16 */
  res = Q2*sqrt(rmax/(2.0*Len*Len*Len)) * exp(-SQR(alpha)*rmax2) / (SQR(alpha)*rmax2);
  
  return res;
}

FLOAT_TYPE compute_error_estimate_k(FLOAT_TYPE alpha, int NP) {
  /* compute the k space part of the error estimate */
  FLOAT_TYPE res;

  /* Kolafa Perram, eq. 31 */
/*   res = Q2 * alpha * pow(PI, -2.0) * pow(kmax, -1.5)  */
/*     * exp(-SQR(PI*kmax/(alpha*L))); */

  /* Petersen 1995, eq. 11 */
  res = 2.0 * Q2 * alpha * Leni * sqrt(1.0/(PI*kmax*NP))
    * exp(-SQR(PI*kmax/(alpha*Len)));

  return res;
}

FLOAT_TYPE Ewald_estimate_error(FLOAT_TYPE alpha, FLOAT_TYPE rmax, int NP) {
  return sqrt(SQR(compute_error_estimate_r(alpha, rmax)) + SQR(compute_error_estimate_k(alpha, NP)));
}

FLOAT_TYPE Ewald_compute_optimal_alpha(FLOAT_TYPE rcut, int NP) {
  /* use bisectional method to get optimal alpha value */
  FLOAT_TYPE alpha_low, f_low;
  FLOAT_TYPE alpha_high, f_high;
  FLOAT_TYPE alpha_guess, f_guess;

  alpha_low = 0.01;
  alpha_high = 10.0;

  f_low = compute_error_estimate_r(alpha_low, rcut) - compute_error_estimate_k(alpha_low, NP);
  f_high = compute_error_estimate_r(alpha_high, rcut) - compute_error_estimate_k(alpha_high, NP);

  if (f_low*f_high > 0.0) {
    printf("Error: Could not init method to find optimal alpha!\n");
    exit(1);
  }

  do {
    alpha_guess = 0.5 *(alpha_low + alpha_high);
    f_guess = compute_error_estimate_r(alpha_guess, rcut) - compute_error_estimate_k(alpha_guess, NP);
    if (f_low*f_guess < 0.0) alpha_high = alpha_guess;
    else alpha_low = alpha_guess;
  } while (fabs(alpha_low-alpha_high) > ALPHA_OPT_PREC);

  return 0.5 *(alpha_low + alpha_high);
}

double Ewald_error_wrapper(double a, int *b, int c, int NP, double e, double alpha_L, double r_cut_iL, double *box_l) {
  return Ewald_estimate_error( alpha_L / box_l[0], r_cut_iL * box_l[0], NP);
}
