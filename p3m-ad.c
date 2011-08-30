
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <string.h>

#include "p3m.h"
#include "charge-assign.h"


const method_t method_p3m_ad = { METHOD_P3M_ad, "P3M with analytic differentiation, not intelaced.", METHOD_FLAG_ad, &Init_ad, &Influence_function_berechnen_ad, &P3M_ad, NULL };

fftw_plan forward_plan;
fftw_plan backward_plan;

FLOAT_TYPE *F_K[3];

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
  
  data_t *d = Init_data( s, p);
  
  forward_plan = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)d->Qmesh, (fftw_complex *)d->Qmesh, FFTW_FORWARD, FFTW_ESTIMATE);

  backward_plan = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)d->Qmesh, (fftw_complex *)d->Qmesh, FFTW_BACKWARD, FFTW_ESTIMATE);
}


void Aliasing_sums_ad(int NX, int NY, int NZ, FLOAT_TYPE alpha,
		      FLOAT_TYPE *Zaehler, FLOAT_TYPE *Nenner1, FLOAT_TYPE *Nenner2)
{
  static int aliasmax = 0; 
  
  FLOAT_TYPE S1,S2,S3;
  FLOAT_TYPE fak1,fak2,zwi;
  int    MX,MY,MZ;
  FLOAT_TYPE NMX,NMY,NMZ;
  FLOAT_TYPE NM2;
  FLOAT_TYPE expo, TE;

  fak1 = 1.0/(FLOAT_TYPE)Mesh;
  fak2 = SQR(PI/(alpha*Len));

  *Zaehler = *Nenner1 = *Nenner2 = 0.0;

  for (MX = -aliasmax; MX <= aliasmax; MX++) {
    NMX = nshift[NX] + Mesh*MX;
    S1   = pow(sinc(fak1*NMX), 2.0*(ip+1)); 
    for (MY = -aliasmax; MY <= aliasmax; MY++) {
      NMY = nshift[NY] + Mesh*MY;
      S2   = S1*pow(sinc(fak1*NMY), 2.0*(ip+1));
      for (MZ = -aliasmax; MZ <= aliasmax; MZ++) {
	NMZ = nshift[NZ] + Mesh*MZ;
	S3   = S2*pow(sinc(fak1*NMZ), 2.0*(ip+1));

	NM2 = SQR(NMX) + SQR(NMY) + SQR(NMZ);
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
  /*
    Berechnet die influence-function, d.h. sowas wie das Produkt aus
    fouriertransformierter Ladungsverschmierung und fouriertransformierter
    Greenschen Funktion.  (-> HOCKNEY/EASTWOOD optimal influence function!)
    alpha  : Ewald-Parameter.
    ip     : Ordnung des charge assigbment schemes.
  */

  int    NX,NY,NZ;
  FLOAT_TYPE Dnx,Dny,Dnz;
  FLOAT_TYPE dMesh,dMeshi;
  FLOAT_TYPE Zaehler=0.0,Nenner1=0.0, Nenner2=0.0;
  FLOAT_TYPE zwi;
  int ind = 0;
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
		G_hat[ind]=0.0;
              else if ((NX%(Mesh/2) == 0) && (NY%(Mesh/2) == 0) && (NZ%(Mesh/2) == 0))
                G_hat[ind]=0.0;
	      else
		{
		  Aliasing_sums_ad(NX,NY,NZ,alpha,&Zaehler,&Nenner1, &Nenner2);
		  zwi = Zaehler / ( Nenner1 * Nenner2 );

		  G_hat[ind] = 2.0 * zwi / PI;
		}
	    }
	}
    }
}


void P3M_ad( system_t *s, parameters_t *p, data_t *d, forces_t *f )
{
  
  /* Zaehlvariablen: */
  int i, j, k, l; 
  /* Schnelles Modulo: */
  int MESHMASKE;
  /* Hilfsvariablen */
  FLOAT_TYPE T1;

  int Mesh = p->mesh;
  
  double dop;

  MESHMASKE = Mesh-1;

  MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;

  memset(d->Qmesh);
  
  /* chargeassignment */
  assign_charge_and_derivatives();
  
  /* Forward Fast Fourier Transform */
  forward_fft();
 
  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<Mesh; k++)
	{
          int c_index = c_ind(i,j,k);

	  T1 = G_hat[r_ind(i,j,k)];
	  d->Qmesh[c_index] *= T1;
	  d->Qmesh[c_index+1] *= T1;
	}

  /* Durchfuehren der Fourier-Rueck-Transformation: */

  backward_fft();

  /* Force assignment */
  assign_forces_ad( Len/(2.0*Mesh*Mesh), f->fields, s->nparticles, d->Qmesh, 0);

  return;
}
