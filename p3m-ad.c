
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <string.h>

#include "p3m.h"
#include "charge-assign.h"

#define CA_DEBUG

/* Pi, weil man's so oft braucht: */
#define PI 3.14159265358979323846264

#define r_ind(A,B,C) ((A)*Mesh*Mesh + (B)*Mesh + (C))
#define c_ind(A,B,C) (2*Mesh*Mesh*(A)+2*Mesh*(B)+2*(C))

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

static FLOAT_TYPE sinc(FLOAT_TYPE d)
{
#define epsi 0.1

#define c2 -0.1666666666667e-0
#define c4  0.8333333333333e-2
#define c6 -0.1984126984127e-3
#define c8  0.2755731922399e-5

  double PId = PI*d, PId2;

  if (fabs(d)>epsi)
    return sin(PId)/PId;
  else {
    PId2 = SQR(PId);
    return 1.0 + PId2*(c2+PId2*(c4+PId2*(c6+PId2*c8)));
  }
  return 1.0;
}


void Init_ad(int Teilchenzahl) {
  int l;
  Qmesh = (FLOAT_TYPE *) realloc(Qmesh, 2*Mesh*Mesh*Mesh*sizeof(FLOAT_TYPE));

  G_hat = (FLOAT_TYPE *) realloc(G_hat, Mesh*Mesh*Mesh*sizeof(FLOAT_TYPE));  

  F_K[0] = Fx_K;
  F_K[1] = Fy_K;
  F_K[2] = Fz_K;

  ca_ind[0] = (int *) realloc(ca_ind[0], 3*Teilchenzahl*sizeof(int));
  cf[0] = (FLOAT_TYPE *) realloc(cf[0], Teilchenzahl * (ip+1) * (ip + 1) * (ip + 1) *sizeof(FLOAT_TYPE)); 

  
  dQdx[0] = (FLOAT_TYPE *) realloc(dQdx[0], Teilchenzahl * (ip+1) * (ip + 1) * (ip + 1) *sizeof(FLOAT_TYPE));
  dQdy[0] = (FLOAT_TYPE *) realloc(dQdy[0], Teilchenzahl * (ip+1) * (ip + 1) * (ip + 1) *sizeof(FLOAT_TYPE));  
  dQdz[0] = (FLOAT_TYPE *) realloc(dQdz[0], Teilchenzahl * (ip+1) * (ip + 1) * (ip + 1) *sizeof(FLOAT_TYPE));

  forward_plan = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)Qmesh, (fftw_complex *)Qmesh, FFTW_FORWARD, FFTW_ESTIMATE);

  backward_plan = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)Qmesh, (fftw_complex *)Qmesh, FFTW_BACKWARD, FFTW_ESTIMATE);
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

void Influence_function_berechnen_ad(FLOAT_TYPE alpha)
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


void P3M_ad(const FLOAT_TYPE alpha, const int Teilchenzahl)
{
  /*
    Berechnet den k-Raum Anteil der Ewald-Routine auf dem Gitter.
    Die Ableitung wird durch analytische Differentiation der charge assignment
    function erreicht, so wie es auch im EPBDLP-paper geschieht.
    alpha : Ewald-Parameter.
  */
  
  /* Zaehlvariablen: */
  int i, j, k, l; 
  /* Variablen fuer FFT: */
  int Nx,Ny,Nz,Lda, Ldb, sx, sy, sz;   
  /* Schnelles Modulo: */
  int MESHMASKE;
  /* Hilfsvariablen */
  FLOAT_TYPE H,Hi,dMesh,MI2;
  /* charge-assignment beschleunigen */
  FLOAT_TYPE T1;
  /* Soweit links vom Referenzpunkt gehts beim Ladungsver-
     teilen los (implementiert ist noch ein Summand Mesh!): */
  FLOAT_TYPE dTeilchenzahli;
  
  int direction;
  double dop;

  Nx = Ny = Nz = Mesh;
  Lda = Ldb = Mesh; 
  sx = sy = sz = 1; 
  MESHMASKE = Mesh-1;
  dMesh = (FLOAT_TYPE)Mesh;
  H = Len/dMesh;
  Hi = 1.0/H;
  MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;
  dTeilchenzahli = 1.0/(FLOAT_TYPE)Teilchenzahl;

  /* Initialisieren von Qmesh */
  bzero(Qmesh, 2*Mesh*Mesh*Mesh*sizeof(FLOAT_TYPE));

  /* chargeassignment */
  for (i=0; i<Teilchenzahl; i++)
    {
      FLOAT_TYPE Ri[3] = {xS[i], yS[i], zS[i]};
      assign_charge_and_derivatives(i, Q[i], Ri, Qmesh, 0);
    }

  /* Forward Fast Fourier Transform */
  forward_fft();
 
  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<Mesh; k++)
	{
          int c_index = c_ind(i,j,k);

	  T1 = G_hat[r_ind(i,j,k)];
	  Qmesh[c_index] *= T1;
	  Qmesh[c_index+1] *= T1;
	}

  /* Durchfuehren der Fourier-Rueck-Transformation: */

  backward_fft();

  /* Force assignment */
  assign_forces_ad( Len/(2.0*Mesh*Mesh), F_K, Teilchenzahl, Qmesh, 0);

  return;
}
