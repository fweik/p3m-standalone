
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <string.h>

#include "p3m.h"

#include "p3m-ad-i.h"

#include "charge-assign.h"

const method_t method_p3m_ad_i = { METHOD_P3M_ad_i, "P3M with analytic differentiation, not intelaced.", 
				   METHOD_FLAG_ad | METHOD_FLAG_interlaced | MEHOD_FLAG_ca | METHOD_FLAG_G_hat, 
				   &Init_ad_i, &Influence_function_berechnen_ad_i, &P3M_ad_i, NULL };

fftw_plan forward_plan;
fftw_plan backward_plan;

static void forward_fft(void);
static void backward_fft(void);

inline void forward_fft(void) {
  fftw_execute(forward_plan);
}

inline void backward_fft(void) {
    fftw_execute(backward_plan);
}

data_t *Init_ad_interlaced( system_t *s, parameters_t *p ) {
  int Mesh = p->mesh;
  data_t *d = Init_ta( &method_p3m_ad_i, s, p );

  forward_plan = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)d->Qmesh, (fftw_complex *)d->Qmesh, FFTW_FORWARD, FFTW_ESTIMATE);

  backward_plan = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)d->Qmesh, (fftw_complex *)d->Qmesh, FFTW_BACKWARD, FFTW_ESTIMATE);
}


void Aliasing_sums_ad_interlaced(int NX, int NY, int NZ, system_t *s, parameters_t *p, data_t *d,
				 FLOAT_TYPE *Zaehler, FLOAT_TYPE *Nenner1, FLOAT_TYPE *Nenner2, FLOAT_TYPE *Nenner3, FLOAT_TYPE *Nenner4)
{
  FLOAT_TYPE S1,S2,S3;
  FLOAT_TYPE fak1,fak2,zwi;
  int    MX,MY,MZ;
  FLOAT_TYPE NMX,NMY,NMZ;
  FLOAT_TYPE NM2;
  FLOAT_TYPE expo, TE;
  
  int Mesh = p->mesh;
  
  fak1 = 1.0/(FLOAT_TYPE)Mesh;
  fak2 = SQR(PI/(p->alpha*s->length));

  *Zaehler = *Nenner1 = *Nenner2 = 0.0;

  for (MX = -P3M_BRILLOUIN; MX <= P3M_BRILLOUIN; MX++) {
    NMX = d->nshift[NX] + Mesh*MX;
    S1   = pow(sinc(fak1*NMX), 2.0*p->cao); 
    for (MY = -P3M_BRILLOUIN; MY <= P3M_BRILLOUIN; MY++) {
      NMY = d->nshift[NY] + Mesh*MY;
      S2   = S1*pow(sinc(fak1*NMY), 2.0*p->cao);
      for (MZ = -P3M_BRILLOUIN; MZ <= P3M_BRILLOUIN; MZ++) {
	NMZ = d->nshift[NZ] + Mesh*MZ;
	S3   = S2*pow(sinc(fak1*NMZ), 2.0*p->cao);

	NM2 = SQR(NMX) + SQR(NMY) + SQR(NMZ);
	*Nenner1 += S3;
	*Nenner2 += S3 * NM2;

        *Nenner3 += ( 1.0 - 2.0*((MX + MY + MZ) % 2)) * S3;
        *Nenner4 += ( 1.0 - 2.0*((MX + MY + MZ) % 2)) * S3 * NM2;

	expo = fak2*NM2;
	TE = (expo < 30) ? exp(-expo) : 0.0;
	zwi  = S3 * TE;
        *Zaehler += zwi;
      }
    }
  }
}

void Influence_function_ad( system_t *s, parameters_t *p, data_t *d )
{
  int    NX,NY,NZ;
  FLOAT_TYPE Dnx,Dny,Dnz;
  FLOAT_TYPE dMesh,dMeshi;
  FLOAT_TYPE Zaehler=0.0,Nenner1=0.0, Nenner2=0.0, Nenner3=0.0, Nenner4=0.0;
  FLOAT_TYPE zwi;
  int ind = 0;
  int Mesh= d->mesh;
  dMesh = (FLOAT_TYPE)Mesh;
  dMeshi= 1.0/dMesh;

  /* bei Zahlen >= Mesh/2 wird noch Mesh abgezogen! */
  for (NX=0; NX<Mesh; NX++)
    {
      for (NY=0; NY<Mesh; NY++)
	{
	  for (NZ=0; NZ<Mesh; NZ++)
	    {
              ind = r_ind(NX,NY,NZ);

	      if ((NX==0) && (NY==0) && (NZ==0))
		d->G_hat[ind]=0.0;
              else if ((NX%(Mesh/2) == 0) && (NY%(Mesh/2) == 0) && (NZ%(Mesh/2) == 0))
                d->G_hat[ind]=0.0;
	      else
		{
		  Aliasing_sums_ad(NX,NY,NZ,p->alpha,&Zaehler,&Nenner1, &Nenner2, &Nenner3, &Nenner4);
		  zwi = Zaehler / ( 0.5 * (Nenner1 * Nenner2 + Nenner3 * Nenner4 ));

		  d->G_hat[ind] = zwi * s->length * s->length / PI;
		}
	    }
	}
    }
}


void P3M_ad_interlaced( system_t *s, parameters_t *p, data_t *d, forces_t *f )
{
  /* Zaehlvariablen: */
  int i, j, k, l, ii; 
  /* Schnelles Modulo: */
  int MESHMASKE;
  /* Hilfsvariablen */
  FLOAT_TYPE H,Hi,dMesh,MI2;
  /* charge-assignment beschleunigen */
  FLOAT_TYPE T1;
  
  int direction;
  double dop;

  int Mesh = p->mesh;

  MESHMASKE = Mesh-1;
  dMesh = (FLOAT_TYPE)Mesh;
  H = s->length/dMesh;
  Hi = 1.0/H;
  MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;

  /* Initialisieren von Qmesh */
  memset(d->Qmesh, 2*Mesh*Mesh*Mesh, sizeof(FLOAT_TYPE));

  /* chargeassignment */
  assign_charges_and_derivatives( s, p, d, 0 );
  assign_charges_and_derivatives( s, p, d, 1 );

  /* Forward Fast Fourier Transform */
  forward_fft();
 
  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<Mesh; k++)
	{
          int c_index = c_ind(i,j,k);

	  T1 = d->G_hat[r_ind(i,j,k)];
	  d->Qmesh[c_index] *= T1;
	  d->Qmesh[c_index+1] *= T1;
	}

  /* Durchfuehren der Fourier-Rueck-Transformation: */

  backward_fft();

  /* Force assignment */
  assign_forces_ad(0.5*H*H*H, s, p, d, f, 0);
  assign_forces_ad(0.5*H*H*H, s, p, d, f, 1);
  return;
}
