
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <string.h>

#include "p3m.h"
#include "charge-assign.h"
#include "p3m-ad.h"

const method_t method_p3m_ad = { METHOD_P3M_ad, "P3M with analytic differentiation, not intelaced.", 
				 METHOD_FLAG_P3M | METHOD_FLAG_ad, 
				 &Init_ad, &Influence_function_berechnen_ad, &P3M_ad, NULL };

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

data_t *Init_ad( system_t *s, parameters_t *p ) {
  int Mesh = p->mesh;
  
  data_t *d = Init_data( &method_p3m_ad, s, p);
  
  forward_plan = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)d->Qmesh, (fftw_complex *)d->Qmesh, FFTW_FORWARD, FFTW_ESTIMATE);

  backward_plan = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)d->Qmesh, (fftw_complex *)d->Qmesh, FFTW_BACKWARD, FFTW_ESTIMATE);

  return d;
}


void Aliasing_sums_ad(int NX, int NY, int NZ, system_t *s, parameters_t *p, data_t *d,
		      FLOAT_TYPE *Zaehler, FLOAT_TYPE *Nenner1, FLOAT_TYPE *Nenner2)
{
  FLOAT_TYPE S1,S2,S3;
  FLOAT_TYPE fak1,fak2,zwi;
  int    MX,MY,MZ;
  FLOAT_TYPE NMX,NMY,NMZ;
  FLOAT_TYPE NM2;
  FLOAT_TYPE expo, TE;
  int Mesh = p->mesh;
  FLOAT_TYPE Len = s->length;
  FLOAT_TYPE Leni = 1.0/Len;

  fak1 = 1.0/(FLOAT_TYPE)Mesh;
  fak2 = SQR(PI/p->alpha);

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

	NM2 = SQR(NMX*Leni) + SQR(NMY*Leni) + SQR(NMZ*Leni);
	*Nenner1 += S3;
	*Nenner2 += S3 * NM2;

	expo = fak2*NM2;
	TE = (expo < 30) ? exp(-expo) : 0.0;
	zwi  = S3 * TE;
        *Zaehler += zwi;
      }
    }
  }
}

void Influence_function_berechnen_ad( system_t *s, parameters_t *p, data_t *d )
{

  int    NX,NY,NZ;
  FLOAT_TYPE dMesh,dMeshi;
  FLOAT_TYPE Zaehler=0.0,Nenner1=0.0, Nenner2=0.0;

  int ind = 0;
  int Mesh = p->mesh;
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
		  Aliasing_sums_ad(NX,NY,NZ,s,p,d,&Zaehler,&Nenner1, &Nenner2);
		  d->G_hat[ind] = Zaehler / ( PI * Nenner1 * Nenner2 );
		}
	    }
	}
    }
}


void P3M_ad( system_t *s, parameters_t *p, data_t *d, forces_t *f )
{
  
  /* Zaehlvariablen: */
  int i, j, k, c_index; 
  /* Hilfsvariablen */
  FLOAT_TYPE T1;
  FLOAT_TYPE Leni = 1.0/s->length;
  int Mesh = p->mesh;
  
  memset(d->Qmesh, 0, 2*Mesh*Mesh*Mesh * sizeof(FLOAT_TYPE));
  
  /* chargeassignment */
  assign_charge_and_derivatives( s, p, d, 0, 1 );
  
  /* Forward Fast Fourier Transform */
  forward_fft();
 
  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<Mesh; k++)
	{
          c_index = c_ind(i,j,k);

	  T1 = d->G_hat[r_ind(i,j,k)];
	  d->Qmesh[c_index] *= T1;
	  d->Qmesh[c_index+1] *= T1;
	}

  /* Durchfuehren der Fourier-Rueck-Transformation: */

  backward_fft();

  /* Force assignment */
  assign_forces_ad( Mesh * Leni * Leni * Leni , s, p, d, f, 0 );

  return;
}
