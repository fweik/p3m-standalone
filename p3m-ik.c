
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


// declaration of the method

const method_t method_p3m_ik = { METHOD_P3M_ik, "P3M with ik differentiation, not intelaced.",
                                 METHOD_FLAG_P3M | METHOD_FLAG_ik,
                                 &Init_ik, &Influence_function_berechnen_ik, &P3M_ik, &Error_ik
                               };

// Forward declaration of local functions

static void forward_fft ( data_t * );
static void backward_fft ( data_t * );
static void p3m_tune_aliasing_sums_ik ( int, int, int,
                                        const system_t *, const parameters_t *,
                                        double *, double * );

static FLOAT_TYPE p3m_k_space_error_ik ( FLOAT_TYPE, const system_t *, const parameters_t * );

inline void forward_fft ( data_t *d ) {
    fftw_execute ( d->forward_plan[0] );
}

inline void backward_fft ( data_t *d ) {
    int i;
    for ( i=0;i<3;i++ )
        fftw_execute ( d->backward_plan[i] );
}

data_t *Init_ik ( system_t *s, parameters_t *p ) {
    int l;
    int mesh = p->mesh;

    data_t *d = Init_data ( &method_p3m_ik, s, p );

    d->forward_plans = 1;
    d->backward_plans = 3;

    d->forward_plan[0] = fftw_plan_dft_3d ( mesh, mesh, mesh, ( fftw_complex * ) d->Qmesh, ( fftw_complex * ) d->Qmesh, FFTW_FORWARD, FFTW_ESTIMATE );

    for ( l=0;l<3;l++ ) {
        d->backward_plan[l] = fftw_plan_dft_3d ( mesh, mesh, mesh, ( fftw_complex * ) ( d->Fmesh->fields[l] ), ( fftw_complex * ) ( d->Fmesh->fields[l] ), FFTW_BACKWARD, FFTW_ESTIMATE );
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

    fak1 = 1.0/ ( FLOAT_TYPE ) Mesh;
    fak2 = SQR ( PI/ ( p->alpha ) );

    Zaehler[0] = Zaehler[1] = Zaehler[2] = *Nenner = 0.0;

    for ( MX = -P3M_BRILLOUIN; MX <= P3M_BRILLOUIN; MX++ ) {
        NMX = d->nshift[NX] + Mesh*MX;
        S1 = pow ( sinc ( fak1*NMX ), 2*p->cao );
        for ( MY = -P3M_BRILLOUIN; MY <= P3M_BRILLOUIN; MY++ ) {
            NMY = d->nshift[NY] + Mesh*MY;
            S2   = S1*pow ( sinc ( fak1*NMY ), 2*p->cao );
            for ( MZ = -P3M_BRILLOUIN; MZ <= P3M_BRILLOUIN; MZ++ ) {
                NMZ = d->nshift[NZ] + Mesh*MZ;
                S3   = S2*pow ( sinc ( fak1*NMZ ), 2*p->cao );

                NM2 = SQR ( NMX*Leni ) + SQR ( NMY*Leni ) + SQR ( NMZ*Leni );
                *Nenner += S3;

                expo = fak2*NM2;
                TE = ( expo < 30.0 ) ? exp ( -expo ) : 0.0;
                zwi  = S3 * TE/NM2;
                Zaehler[0] += NMX*zwi*Leni;
                Zaehler[1] += NMY*zwi*Leni;
                Zaehler[2] += NMZ*zwi*Leni;
            }
        }
    }
}

void Influence_function_berechnen_ik ( system_t *s, parameters_t *p, data_t *d ) {

  int    NX,NY,NZ;
  FLOAT_TYPE Dnx,Dny,Dnz;
  FLOAT_TYPE dMesh,dMeshi;
  FLOAT_TYPE Zaehler[3]={0.0,0.0,0.0},Nenner=0.0;
  FLOAT_TYPE zwi;
  int ind = 0;
  int Mesh = p->mesh;
  FLOAT_TYPE Leni = 1.0/s->length;
  dMesh = ( FLOAT_TYPE ) Mesh;
  dMeshi= 1.0/dMesh;

  if(Mesh > 512) 
    exit(-23);

  for ( NX=0; NX<Mesh; NX++ ) {
    for ( NY=0; NY<Mesh; NY++ ) {
      for ( NZ=0; NZ<Mesh; NZ++ ) {
	ind = r_ind ( NX, NY, NZ );
	  
	if ( ( NX==0 ) && ( NY==0 ) && ( NZ==0 ) )
	  d->G_hat[ind]=0.0;
	else if ( ( NX% ( Mesh/2 ) == 0 ) && ( NY% ( Mesh/2 ) == 0 ) && ( NZ% ( Mesh/2 ) == 0 ) )
	  d->G_hat[ind]=0.0;
	else {
	  Aliasing_sums_ik ( s, p, d, NX,NY,NZ,Zaehler,&Nenner );
		  
	  Dnx = d->Dn[NX];
	  Dny = d->Dn[NY];
	  Dnz = d->Dn[NZ];
	    
	  zwi  = Dnx*Zaehler[0]*Leni + Dny*Zaehler[1]*Leni + Dnz*Zaehler[2]*Leni;
	  zwi /= ( ( SQR ( Dnx*Leni ) + SQR ( Dny*Leni ) + SQR ( Dnz*Leni ) ) * SQR ( Nenner ) );
	  d->G_hat[ind] = 2.0 * zwi / PI;
	}
      }
    }
  }
}


/* Calculates k-space part of the force, using ik-differentiation.
 */

void P3M_ik ( system_t *s, parameters_t *p, data_t *d, forces_t *f ) {
    /* Zaehlvariablen: */
    int i, j, k, l;
    /* Hilfsvariablen */
    FLOAT_TYPE H,Hi,dMesh,MI2;
    /* charge-assignment beschleunigen */
    FLOAT_TYPE T1;
    /* Soweit links vom Referenzpunkt gehts beim Ladungsver-
       teilen los (implementiert ist noch ein Summand Mesh!): */
    FLOAT_TYPE dTeilchenzahli;
    FLOAT_TYPE Leni = 1.0/s->length;
    double dop;
    int Mesh = p->mesh;
    int c_index;

    dMesh = ( FLOAT_TYPE ) Mesh;
    H = s->length/dMesh;
    Hi = 1.0/H;
    MI2 = 2.0* ( FLOAT_TYPE ) MaxInterpol;
    dTeilchenzahli = 1.0/ ( FLOAT_TYPE ) s->nparticles;

    /* Initialisieren von Qmesh */
    memset ( d->Qmesh, 0, 2*Mesh*Mesh*Mesh*sizeof ( FLOAT_TYPE ) );

    /* chargeassignment */
    assign_charge ( s, p, d, 0 );

    /* Forward Fast Fourier Transform */
    forward_fft(d);

    for ( i=0; i<Mesh; i++ )
        for ( j=0; j<Mesh; j++ )
            for ( k=0; k<Mesh; k++ ) {
                c_index = c_ind ( i,j,k );

                T1 = d->G_hat[r_ind ( i,j,k ) ];
                d->Qmesh[c_index] *= T1;
                d->Qmesh[c_index+1] *= T1;
 
                for ( l=0;l<3;l++ ) {
                    switch ( l ) {
                    case 0:
                        dop = d->Dn[i];
                        break;
                    case 1:
                        dop = d->Dn[j];
                        break;
                    case 2:
                        dop = d->Dn[k];
                        break;
                    }
                    d->Fmesh->fields[l][c_index]   =  -2.0*PI*Leni*dop*d->Qmesh[c_index+1];
                    d->Fmesh->fields[l][c_index+1] =   2.0*PI*Leni*dop*d->Qmesh[c_index];
                }
            }

    /* Durchfuehren der Fourier-Rueck-Transformation: */
    backward_fft(d);

    /* Force assignment */
    assign_forces ( 1.0/ ( 2.0*s->length*s->length*s->length ),s,p,d,f,0 );
    return;
}

// Functions for error estimate.

FLOAT_TYPE Error_ik( system_t *s, parameters_t *p) {
  return sqrt( SQR( Realspace_error( s, p ) ) + SQR( p3m_k_space_error_ik( 1.0, s, p) ) );
}

FLOAT_TYPE p3m_k_space_error_ik ( FLOAT_TYPE prefac, const system_t *s, const parameters_t *p ) {
    int  nx, ny, nz;
    FLOAT_TYPE he_q = 0.0;
    FLOAT_TYPE alias1, alias2, n2, cs;
    FLOAT_TYPE ctan_x, ctan_y;
    int mesh = p->mesh;
    FLOAT_TYPE meshi = 1.0/(FLOAT_TYPE)(p->mesh);

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
    return 2.0*s->q2*sqrt ( he_q/ ( FLOAT_TYPE ) s->nparticles ) / ( SQR ( s->length ) );
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
                ex = exp ( -factor1*nm2 );
                ex2 = SQR ( ex );

                U2 = pow ( sinc ( fnmx ) *sinc ( fnmy ) *sinc ( fnmz ), 2.0*p->cao );

                *alias1 += ex2 / nm2;
                *alias2 += U2 * ex * ( nx*nmx + ny*nmy + nz*nmz ) / nm2;
            }
        }
    }
}
