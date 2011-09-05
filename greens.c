
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <string.h>

#include "types.h"
#include "common.h"
#include "p3m-common.h"
#include "charge-assign.h"
#include "greens.h"

fftw_plan forward_plan;
fftw_plan backward_plan[3];

// declaration of the method

const method_t method_greens_ik = { METHOD_GREENS_ik, "Continuum Greens function with ik differentiation, not intelaced.",
                                 METHOD_FLAG_ik | METHOD_FLAG_Qmesh | METHOD_FLAG_G_hat | METHOD_FLAG_nshift | METHOD_FLAG_ca,
                                 &Init_greens_ik, &Greens_function, &Greens_kspace_ik, NULL
                               };

// Forward declaration of local functions

static void forward_fft ( void );
static void backward_fft ( void );

inline void forward_fft ( void ) {
    fftw_execute ( forward_plan );
}

inline void backward_fft ( void ) {
    int i;
    for ( i=0;i<3;i++ )
        fftw_execute ( backward_plan[i] );
}

data_t *Init_greens_ik ( system_t *s, parameters_t *p ) {
    int l;
    int mesh = p->mesh;

    data_t *d = Init_data ( &method_greens_ik, s, p );

    forward_plan = fftw_plan_dft_3d ( mesh, mesh, mesh, ( fftw_complex * ) d->Qmesh, ( fftw_complex * ) d->Qmesh, FFTW_FORWARD, FFTW_ESTIMATE );

    for ( l=0;l<3;l++ ) {
        backward_plan[l] = fftw_plan_dft_3d ( mesh, mesh, mesh, ( fftw_complex * ) ( d->Fmesh->fields[l] ), ( fftw_complex * ) ( d->Fmesh->fields[l] ), FFTW_BACKWARD, FFTW_ESTIMATE );
    }
    return d;
}

void Greens_function ( system_t *s, parameters_t *p, data_t *d ) {

    int    NX,NY,NZ;
    FLOAT_TYPE Dnx,Dny,Dnz;
    FLOAT_TYPE dMesh,dMeshi;
    FLOAT_TYPE NM2, expo, TE;
    int ind = 0;
    int Mesh = p->mesh;
    FLOAT_TYPE Leni = 1.0/s->length;
    dMesh = ( FLOAT_TYPE ) Mesh;
    dMeshi= 1.0/dMesh;
    
    for ( NX=0; NX<Mesh; NX++ ) {
      for ( NY=0; NY<Mesh; NY++ ) {
	for ( NZ=0; NZ<Mesh; NZ++ ) {
	  ind = r_ind ( NX,NY,NZ );
	  
	  if ( ( ( NX==0 ) && ( NY==0 ) && ( NZ==0 ) )
            || ( ( NX% ( Mesh/2 ) == 0 ) && ( NY% ( Mesh/2 ) == 0 ) && ( NZ% ( Mesh/2 ) == 0 ) ) )
	    d->G_hat[ind]=0.0;
	  else {
		  
	    Dnx = d->Dn[NX];
	    Dny = d->Dn[NY];
	    Dnz = d->Dn[NZ];
	    
            NM2 = SQR ( Dnx*Leni ) + SQR ( Dny*Leni ) + SQR ( Dnz*Leni );
            
            expo = NM2 * SQR ( PI / p->alpha );
            TE = ( expo < 30.0 ) ? exp ( -expo ) : 0.0;
            
	    d->G_hat[ind] = TE / ( SQR(2.0*PI)*NM2 );
	  }
	}
      }
    }
}


/* Calculates k-space part of the force, using ik-differentiation.
 */

void Greens_kspace_ik ( system_t *s, parameters_t *p, data_t *d, forces_t *f ) {
    /* Zaehlvariablen: */
    int i, j, k, l;
    /* Schnelles Modulo: */
    int MESHMASKE;
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

    MESHMASKE = Mesh-1;
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
    forward_fft();

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
    backward_fft();

    /* Force assignment */
    assign_forces ( 4.0*PI / (s->length * s->length * s->length),  s,p,d,f,0 );
    return;
}
