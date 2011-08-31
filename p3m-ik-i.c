
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

#include "p3m.h"
#include "common.h"
#include "p3m-common.h"
#include "charge-assign.h"
#include "p3m-ik-i.h"

// declaration of the method

const method_t method_p3m_ik_i = { METHOD_P3M_ik_i, "P3M with ik differentiation, interlaced.",
                                   METHOD_FLAG_ik | METHOD_FLAG_Qmesh | METHOD_FLAG_G_hat | METHOD_FLAG_nshift | METHOD_FLAG_ca | METHOD_FLAG_interlaced,
                                   &Init_ik_i, &Influence_ik_i, &P3M_ik_i, NULL
                                 };

fftw_plan forward_plan;
fftw_plan backward_plan[3];

static void forward_fft(void);
static void backward_fft(void);

void forward_fft(void) {
    fftw_execute(forward_plan);
}

void backward_fft(void) {
    int i;
    for (i=0;i<3;i++)
        fftw_execute(backward_plan[i]);
}

data_t *Init_ik_i( system_t *s, parameters_t *p ) {
    int l;
    int Mesh = p->mesh;

    data_t *d = Init_data( &method_p3m_ik_i, s, p );

    forward_plan = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)d->Qmesh, (fftw_complex *)d->Qmesh, FFTW_FORWARD, FFTW_ESTIMATE);
    for (l=0;l<3;l++) {
        backward_plan[l] = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)d->Fmesh->fields[l], (fftw_complex *)d->Fmesh->fields[l], FFTW_BACKWARD, FFTW_ESTIMATE);
    }
    return d;
}

void Aliasing_sums_ik_i( system_t *s, parameters_t *p, data_t *d, int NX, int NY, int NZ,
                         FLOAT_TYPE *Zaehler, FLOAT_TYPE *Nenner1, FLOAT_TYPE *Nenner2)  {
  FLOAT_TYPE S1,S2,S3;
  FLOAT_TYPE fak1,fak2,zwi;
  int    MX,MY,MZ;
  FLOAT_TYPE NMX,NMY,NMZ;
  FLOAT_TYPE NM2;
  
  FLOAT_TYPE expo, TE;
  int Mesh = d->mesh;
    
  
  
  
  FLOAT_TYPE Leni = 1.0/s->length;
  fak1 = 1.0/(FLOAT_TYPE)Mesh;
  fak2 = SQR(PI/p->alpha);
  
  Zaehler[0] = Zaehler[1] = Zaehler[2] = *Nenner1 = *Nenner2 = 0.0;
  
  for ( MX = -P3M_BRILLOUIN; MX <= P3M_BRILLOUIN; MX++ ) {
    NMX = d->nshift[NX] + Mesh*MX;
    S1  = pow( sinc(fak1*NMX), 2*p->cao );
    for ( MY = -P3M_BRILLOUIN; MY <= P3M_BRILLOUIN; MY++ ) {
      NMY = d->nshift[NY] + Mesh*MY;
      S2   = S1 * pow( sinc(fak1*NMY), 2*p->cao );
      for ( MZ = -P3M_BRILLOUIN; MZ <= P3M_BRILLOUIN; MZ++ ) {
	NMZ = d->nshift[NZ] + Mesh*MZ;
	S3  = S2*pow( sinc(fak1*NMZ), 2*p->cao );
	
	NM2 = SQR ( NMX*Leni ) + SQR ( NMY*Leni ) + SQR ( NMZ*Leni );
	*Nenner1 += S3;
	
	expo = fak2*NM2;
	TE = ( expo < 30.0 ) ? exp ( -expo ) : 0.0;
	zwi  = S3 * TE/NM2;
	Zaehler[0] += NMX*zwi*Leni;
	Zaehler[1] += NMY*zwi*Leni;
	Zaehler[2] += NMZ*zwi*Leni;
	
	if (((MX+MY+MZ)%2)==0)                                       //even term
	  *Nenner2 += S3;
	else                                                //odd term: minus sign!
	  *Nenner2 -= S3;
      }
    }
  }
}



void Influence_ik_i( system_t *s, parameters_t *p, data_t *d )
{
    int    NX,NY,NZ;
    FLOAT_TYPE Dnx,Dny,Dnz;
    FLOAT_TYPE dMesh,dMeshi;
    FLOAT_TYPE Zaehler[3]={0.0,0.0,0.0},Nenner1=0.0, Nenner2=0.0;
    FLOAT_TYPE zwi;

    int ind = 0;
    int Mesh = p->mesh;
    FLOAT_TYPE Leni = 1.0/s->length;

    dMesh = (FLOAT_TYPE)Mesh;
    dMeshi= 1.0/dMesh;

    for (NX=0; NX<Mesh; NX++)
    {
        for (NY=0; NY<Mesh; NY++)
        {
            for (NZ=0; NZ<Mesh; NZ++)
            {
                ind = r_ind( NX, NY, NZ );

                if ((NX==0) && (NY==0) && (NZ==0))
                    d->G_hat[ind]=0.0;
                else if ((NX%(Mesh/2) == 0) && (NY%(Mesh/2) == 0) && (NZ%(Mesh/2) == 0))
                    d->G_hat[ind]=0.0;
                else
                {
                    Aliasing_sums_ik_i( s, p, d, NX, NY, NZ, Zaehler, &Nenner1,  &Nenner2);

                    Dnx = d->Dn[NX];
                    Dny = d->Dn[NY];
                    Dnz = d->Dn[NZ];

                    zwi  = Dnx*Zaehler[0]*Leni + Dny*Zaehler[1]*Leni + Dnz*Zaehler[2]*Leni;
                    zwi /= ( SQR(Dnx*Leni) + SQR(Dny*Leni) + SQR(Dnz*Leni) );
                    zwi /= 0.5*(SQR(Nenner1) + SQR(Nenner2));

                    d->G_hat[ind] = 2.0 * zwi / PI;
                }
            }
        }
    }
}


void P3M_ik_i( system_t *s, parameters_t *p, data_t *d, forces_t *f )
{
    /* Zaehlvariablen: */
    int i, j, k, l;
    /* Schnelles Modulo: */

    FLOAT_TYPE T1;

    FLOAT_TYPE Mesh = p->mesh;
    FLOAT_TYPE Leni = 1.0/s->length;

    FLOAT_TYPE dop;

    int c_index;

    /* Initialisieren von Qmesh */
    memset ( d->Qmesh, 0, 2*Mesh*Mesh*Mesh*sizeof ( FLOAT_TYPE ) );

    /* chargeassignment */
    assign_charge( s, p, d, 0 );
    assign_charge( s, p, d, 1 );

    /* Durchfuehren der Fourier-Hin-Transformationen: */
    forward_fft();

    for (i=0; i<Mesh; i++)
        for (j=0; j<Mesh; j++)
            for (k=0; k<Mesh; k++)
            {
                c_index = c_ind(i,j,k);

                T1 = d->G_hat[r_ind(i,j,k)];
                d->Qmesh[c_index] *= T1;
                d->Qmesh[c_index+1] *= T1;

                for (l=0;l<3;l++) {
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
                    d->Fmesh->fields[l][c_index]   = -2.0*PI*Leni*dop*d->Qmesh[c_index+1];
                    d->Fmesh->fields[l][c_index+1] =  2.0*PI*Leni*dop*d->Qmesh[c_index];
                }

            }

    /* Durchfuehren der Fourier-Rueck-Transformation: */
    backward_fft();

    /* force assignment */
    assign_forces ( 1.0 / ( 2.0*s->length*s->length*s->length ), s, p, d, f, 0 );
    assign_forces ( 1.0 / ( 2.0*s->length*s->length*s->length ), s, p, d, f, 1 );

    return;
}
