
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <string.h>

// General typ definitions
#include "types.h"
#include "common.h"
#include "p3m-common.h"
// Charge assignment
#include "charge-assign.h"

// For realpart error
#include "realpart.h"

#include "p3m-ik.h"

#ifdef __detailed_timings
#include <mpi.h>
#endif


// declaration of the method

const method_t method_p3m_ik = { METHOD_P3M_ik, "P3M with ik differentiation, not intelaced.", "p3m-ik",
                                 METHOD_FLAG_P3M | METHOD_FLAG_ik,
                                 &Init_ik, &Influence_function_berechnen_ik, &P3M_ik, &Error_ik, &Error_ik_k,
                               };

// Forward declaration of local functions

static void forward_fft ( data_t * );
static void backward_fft ( data_t * );
static void p3m_tune_aliasing_sums_ik ( int, int, int,
                                        const system_t *, const parameters_t *,
                                        FLOAT_TYPE *, FLOAT_TYPE * );

FLOAT_TYPE p3m_k_space_error_ik ( FLOAT_TYPE prefac, const system_t *s, const parameters_t *p );

inline void forward_fft ( data_t *d ) {
    FFTW_EXECUTE ( d->forward_plan[0] );
}

inline void backward_fft ( data_t *d ) {
    int i;
    for ( i=0;i<3;i++ )
        FFTW_EXECUTE ( d->backward_plan[i] );
}

FLOAT_TYPE Error_ik_k( system_t *s, parameters_t *p ) {
  return p3m_k_space_error_ik( 1.0, s, p );
}

data_t *Init_ik ( system_t *s, parameters_t *p ) {
    int l;
    int mesh = p->mesh;

    data_t *d = Init_data ( &method_p3m_ik, s, p );

    d->forward_plans = 1;
    d->backward_plans = 3;

    d->forward_plan[0] = FFTW_PLAN_DFT_3D ( mesh, mesh, mesh, ( FFTW_COMPLEX * ) d->Qmesh, ( FFTW_COMPLEX * ) d->Qmesh, FFTW_FORWARD, FFTW_PATIENT );

    for ( l=0;l<3;l++ ) {
        d->backward_plan[l] = FFTW_PLAN_DFT_3D ( mesh, mesh, mesh, ( FFTW_COMPLEX * ) ( d->Fmesh->fields[l] ), ( FFTW_COMPLEX * ) ( d->Fmesh->fields[l] ), FFTW_BACKWARD, FFTW_PATIENT );
    }
    return d;
}


void Aliasing_sums_ik ( system_t *s, parameters_t *p, data_t *d, int NX, int NY, int NZ,
                        FLOAT_TYPE *Zaehler, FLOAT_TYPE *Nenner ) {
    FLOAT_TYPE S1,S2,S3;
    FLOAT_TYPE fak1,fak2,zwi;
    int    MX,MY,MZ;
    FLOAT_TYPE NMX,NMY,NMZ;
    FLOAT_TYPE NM2;
    FLOAT_TYPE expo, TE;
    int Mesh = d->mesh;
    FLOAT_TYPE Leni = 1.0/s->length;
    FLOAT_TYPE (*U_hat)(int, FLOAT_TYPE) = d->inter->U_hat;

    fak1 = 1.0/ ( FLOAT_TYPE ) Mesh;
    fak2 = SQR ( PI/ ( p->alpha ) );

    Zaehler[0] = Zaehler[1] = Zaehler[2] = *Nenner = 0.0;

    for ( MX = -P3M_BRILLOUIN; MX <= P3M_BRILLOUIN; MX++ ) {
        NMX = d->nshift[NX] + Mesh*MX;
        S1 = my_power ( U_hat ( p->cao, fak1*NMX ), 2 );
        for ( MY = -P3M_BRILLOUIN; MY <= P3M_BRILLOUIN; MY++ ) {
            NMY = d->nshift[NY] + Mesh*MY;
            S2   = S1*my_power ( U_hat ( p->cao, fak1*NMY ), 2 );
            for ( MZ = -P3M_BRILLOUIN; MZ <= P3M_BRILLOUIN; MZ++ ) {
                NMZ = d->nshift[NZ] + Mesh*MZ;
                S3   = S2*my_power ( U_hat ( p->cao, fak1*NMZ ), 2 );

                NM2 = SQR ( NMX*Leni ) + SQR ( NMY*Leni ) + SQR ( NMZ*Leni );
                *Nenner += S3;

                expo = fak2*NM2;
                TE = EXP ( -expo );
                zwi  = S3 * TE/NM2;
                Zaehler[0] += NMX*zwi*Leni;
                Zaehler[1] += NMY*zwi*Leni;
                Zaehler[2] += NMZ*zwi*Leni;
            }
        }
    }
}

/* Calculate influence function */
void Influence_function_berechnen_ik ( system_t *s, parameters_t *p, data_t *d ) {

  int    NX,NY,NZ;
  FLOAT_TYPE Dnx,Dny,Dnz;
  FLOAT_TYPE Zaehler[3]={0.0,0.0,0.0},Nenner=0.0;
  FLOAT_TYPE zwi;
  int ind = 0;
  int Mesh = p->mesh;
  FLOAT_TYPE Leni = 1.0/s->length;
  for ( NX=0; NX<Mesh; NX++ ) {
    Dnx = d->Dn[NX];
    for ( NY=0; NY<Mesh; NY++ ) {
      Dny = d->Dn[NY];
#ifdef _OPENMP
#pragma omp parallel for private(NZ, ind, Zaehler, Nenner, Dnz, zwi)
#endif
      for ( NZ=0; NZ<Mesh; NZ++ ) {
	ind = r_ind ( NX, NY, NZ );
	  
	if ( ( NX==0 ) && ( NY==0 ) && ( NZ==0 ) )
	  d->G_hat[ind]=0.0;
	else if ( ( NX% ( Mesh/2 ) == 0 ) && ( NY% ( Mesh/2 ) == 0 ) && ( NZ% ( Mesh/2 ) == 0 ) )
	  d->G_hat[ind]=0.0;
	else {
	  Aliasing_sums_ik ( s, p, d, NX, NY, NZ, Zaehler, &Nenner );
		  
	  Dnz = d->Dn[NZ];
	    
	  zwi  = Dnx*Zaehler[0]*Leni + Dny*Zaehler[1]*Leni + Dnz*Zaehler[2]*Leni;
	  zwi /= ( ( SQR ( Dnx*Leni ) + SQR ( Dny*Leni ) + SQR ( Dnz*Leni ) ) * SQR ( Nenner ) );
	  d->G_hat[ind] = 2.0 * zwi / PI;
	}
      }
    }
  }
#ifdef _OPENMP
#pragma omp barrier
#endif
}


/* Calculates k-space part of the force, using ik-differentiation.
 */

void P3M_ik ( system_t *s, parameters_t *p, data_t *d, forces_t *f ) {
    /* Loop counters */
    int i, j, k;
    /* helper variables */
    FLOAT_TYPE T1;
    FLOAT_TYPE dop;

    // One over boxlength
    FLOAT_TYPE Leni = 1.0/s->length;

    int Mesh = p->mesh;
    int c_index;

    /* Setting charge mesh to zero */
    memset ( d->Qmesh, 0, 2*Mesh*Mesh*Mesh*sizeof ( FLOAT_TYPE ) );

  #ifdef __detailed_timings
  timer = MPI_Wtime();
  #endif


    /* chargeassignment */
    assign_charge ( s, p, d, 0 );

  #ifdef __detailed_timings
  timer = MPI_Wtime() - timer;
  t_charge_assignment[0] = timer;
  timer = MPI_Wtime();
  #endif

    /* Forward Fast Fourier Transform */
    forward_fft(d);

  #ifdef __detailed_timings
    timer = MPI_Wtime() - timer;
    t_fft[0] = timer;
   timer = MPI_Wtime();
  #endif

   double q_r, q_i;


    /* Convolution */
    for ( i=0; i<Mesh; i++ )
        for ( j=0; j<Mesh; j++ )
            for ( k=0; k<Mesh; k++ ) {
                c_index = c_ind ( i,j,k );

                T1 = d->G_hat[r_ind ( i,j,k ) ];
		q_r = -2.0*PI*Leni*d->Qmesh[c_index] *T1;
		q_i = 2.0*PI*Leni*d->Qmesh[c_index+1] *T1;

		dop = d->Dn[i];	 
		d->Fmesh->fields[0][c_index]   =  dop*q_r;
		d->Fmesh->fields[0][c_index+1] =  dop*q_i;

		dop = d->Dn[j];
		d->Fmesh->fields[1][c_index]   =  dop*q_r;
		d->Fmesh->fields[1][c_index+1] =  dop*q_i;

		dop = d->Dn[k];
		d->Fmesh->fields[2][c_index]   =  dop*q_r;
		d->Fmesh->fields[2][c_index+1] =  dop*q_i;
 
            }

  #ifdef __detailed_timings
    timer = MPI_Wtime() - timer;
    t_convolution[0] = timer;
   timer = MPI_Wtime();
  #endif

    /* Backward Fast Fourier Transformation */
    backward_fft(d);

  #ifdef __detailed_timings
    timer = MPI_Wtime() - timer;
    t_fft[0] += timer;
   timer = MPI_Wtime();
  #endif

    /* Force assignment */
    assign_forces ( 1.0/ ( 2.0*s->length*s->length*s->length ),s,p,d,f,0 );

  #ifdef __detailed_timings
    timer = MPI_Wtime() - timer;
    t_force_assignment[0] = timer;
  #endif
}

// Functions for error estimate.

FLOAT_TYPE Error_ik( system_t *s, parameters_t *p) {
  FLOAT_TYPE real;
  FLOAT_TYPE recp;

  real = Realspace_error( s, p );
  recp = p3m_k_space_error_ik( 1.0, s, p);

  //  printf("p3m ik error for mesh %d rcut %lf cao %d alpha %lf : real %e recp %e\n", p->mesh, p->rcut, p->cao, p->alpha, real, recp);
  //printf("system size %d box %lf\n", s->nparticles, s->length);


  return SQRT( SQR( real ) + SQR( recp ) );
}

FLOAT_TYPE p3m_k_space_error_ik ( FLOAT_TYPE prefac, const system_t *s, const parameters_t *p ) {
  // Mesh loop counters 
    int  nx, ny, nz;
    // The pair force error Q
    FLOAT_TYPE he_q = 0.0;
    // Helper variables
    FLOAT_TYPE alias1, alias2, n2, cs;
    FLOAT_TYPE ctan_x, ctan_y;
    int mesh = p->mesh;
    FLOAT_TYPE meshi = 1.0/(FLOAT_TYPE)(p->mesh);

#ifdef _OPENMP
#pragma omp parallel for private(ctan_x, ctan_y, n2, cs, alias1, alias2, ny, nz) reduction( + : he_q )
#endif
    for ( nx=-mesh/2; nx<mesh/2; nx++ ) {
        ctan_x = analytic_cotangent_sum ( nx, meshi,p->cao );
        for ( ny=-mesh/2; ny<mesh/2; ny++ ) {
            ctan_y = ctan_x * analytic_cotangent_sum ( ny, meshi, p->cao );
            for ( nz=-mesh/2; nz<mesh/2; nz++ ) {
                if ( ( nx!=0 ) || ( ny!=0 ) || ( nz!=0 ) ) {
                    n2 = SQR ( nx ) + SQR ( ny ) + SQR ( nz );
                    cs = analytic_cotangent_sum ( nz, meshi ,p->cao ) *ctan_y;
                    p3m_tune_aliasing_sums_ik ( nx,ny,nz, s, p, &alias1,&alias2 );
                    he_q += ( alias1  -  SQR ( alias2/cs ) / n2 );
                }
            }
        }
    }
#ifdef _OPENMP
#pragma omp barrier
#endif
    he_q = fabs(he_q);
    return 2.0*s->q2*SQRT ( he_q/ ( FLOAT_TYPE ) s->nparticles ) / ( SQR ( s->length ) );
}


void p3m_tune_aliasing_sums_ik ( int nx, int ny, int nz,
                                 const system_t *s, const parameters_t *p,
                                 FLOAT_TYPE *alias1, FLOAT_TYPE *alias2 ) {

    int    mx,my,mz;
    FLOAT_TYPE nmx,nmy,nmz;
    FLOAT_TYPE fnmx,fnmy,fnmz;

    FLOAT_TYPE ex,ex2,nm2,U2,factor1;
    int mesh = p->mesh;
    FLOAT_TYPE meshi = 1.0/(FLOAT_TYPE)(p->mesh);
  
    factor1 = SQR ( PI / ( p->alpha * s->length ) );

    *alias1 = *alias2 = 0.0;
    for ( mx=-P3M_BRILLOUIN_TUNING; mx<=P3M_BRILLOUIN_TUNING; mx++ ) {
        fnmx = meshi * ( nmx = nx + mx*mesh );
        for ( my=-P3M_BRILLOUIN_TUNING; my<=P3M_BRILLOUIN_TUNING; my++ ) {
            fnmy = meshi * ( nmy = ny + my*mesh );
            for ( mz=-P3M_BRILLOUIN_TUNING; mz<=P3M_BRILLOUIN_TUNING; mz++ ) {
                fnmz = meshi * ( nmz = nz + mz*mesh );

                nm2 = SQR ( nmx ) + SQR ( nmy ) + SQR ( nmz );
                ex = EXP ( -factor1*nm2 );
                ex2 = SQR ( ex );

                U2 = my_power ( sinc ( fnmx ) *sinc ( fnmy ) *sinc ( fnmz ), 2*p->cao );

                *alias1 += ex2 / nm2;
                *alias2 += U2 * ex * ( nx*nmx + ny*nmy + nz*nmz ) / nm2;
            }
        }
    }
}
