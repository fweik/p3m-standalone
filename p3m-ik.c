
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
fftw_plan backward_plan[3];

FLOAT_TYPE *Fmesh[3];
FLOAT_TYPE *F_K[3];

static void forward_fft(void);
static void backward_fft(void);

inline void forward_fft(void) {
  fftw_execute(forward_plan);
}

inline void backward_fft(void) {
  int i;
  for(i=0;i<3;i++)
    fftw_execute(backward_plan[i]);
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


void Init_ik(int Teilchenzahl) {
  int l;
  Qmesh = (FLOAT_TYPE *) realloc(Qmesh, 2*Mesh*Mesh*Mesh*sizeof(FLOAT_TYPE));

  G_hat = (FLOAT_TYPE *) realloc(G_hat, Mesh*Mesh*Mesh*sizeof(G_hat));  

  F_K[0] = Fx_K;
  F_K[1] = Fy_K;
  F_K[2] = Fz_K;

  ca_ind[0] = (int *) realloc(ca_ind[0], 3*Teilchenzahl*sizeof(int));
  cf[0] = (FLOAT_TYPE *) realloc(cf[0], Teilchenzahl * (ip+1) * (ip + 1) * (ip + 1) *sizeof(FLOAT_TYPE)); 

  forward_plan = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)Qmesh, (fftw_complex *)Qmesh, FFTW_FORWARD, FFTW_ESTIMATE);

  for(l=0;l<3;l++){
    Fmesh[l] = (FLOAT_TYPE *)realloc(Fmesh[l], 2*Mesh*Mesh*Mesh*sizeof(FLOAT_TYPE));  
    backward_plan[l] = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)Fmesh[l], (fftw_complex *)Fmesh[l], FFTW_BACKWARD, FFTW_ESTIMATE);
  }
  
}


void Aliasing_sums_ik(int NX, int NY, int NZ, FLOAT_TYPE alpha,
				  FLOAT_TYPE *Zaehler, FLOAT_TYPE *Nenner)
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

  Zaehler[0] = Zaehler[1] = Zaehler[2] = *Nenner = 0.0;

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
	*Nenner += S3;
	
	expo = fak2*NM2;
	TE = (expo < 30) ? exp(-expo) : 0.0;
	zwi  = S3 * TE/NM2;
	Zaehler[0] += NMX*zwi;
	Zaehler[1] += NMY*zwi;
	Zaehler[2] += NMZ*zwi;
      }
    }
  }
}

void Influence_function_berechnen_ik(FLOAT_TYPE alpha)
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
  FLOAT_TYPE Zaehler[3]={0.0,0.0,0.0},Nenner=0.0;
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
		  Aliasing_sums_ik(NX,NY,NZ,alpha,Zaehler,&Nenner);
                  
		  Dnx = Dn[NX];  
		  Dny = Dn[NY];  
		  Dnz = Dn[NZ];
		  
		  zwi  = Dnx*Zaehler[0] + Dny*Zaehler[1] + Dnz*Zaehler[2];
		  zwi /= ( (SQR(Dnx) + SQR(Dny) + SQR(Dnz)) * SQR(Nenner) );
		  G_hat[ind] = Mesh*Mesh*Mesh*2.0 * zwi * Leni * Leni * Leni;
		}
	    }
	}
    }
}


void P3M_ik(const FLOAT_TYPE alpha, const int Teilchenzahl)
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
      assign_charge(i, Q[i], Ri, Qmesh, 0);
    }

#ifdef CA_DEBUG
  for(i=0;i<Mesh;i++)
    for(j=0;j<Mesh;j++)
      for(k=0;k<Mesh;k++) 
	printf("CA_DEBUG: %lf\n", Qmesh[c_ind(i,j,k)]);
#endif

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

          for(l=0;l<3;l++) {
            switch(l) {
	      case 0:
                dop = Dn[i];
                break;
	      case 1:
                dop = Dn[j];
                break;
	      case 2:
                dop = Dn[k];
                break;
	    }
	    Fmesh[l][c_index]   = -dop*Qmesh[c_index+1];
	    Fmesh[l][c_index+1] =  dop*Qmesh[c_index];
          }
	}

  /* Durchfuehren der Fourier-Rueck-Transformation: */

  backward_fft();
  for(l=0;l<3;l++)
    for(i=0;i<Mesh;i++)
      for(j=0;j<Mesh;j++)
	for(k=0;k<Mesh;k++)
          if(Fmesh[l][c_ind(i,j,k)+1] > 1e-5) {
            printf("WARNING! UNREAL FORCE AT (%d %d %d %d): %lf", i,j,k,l,Fmesh[l][c_ind(i,j,k)+1]);
            exit(128);
	  }

  /* Force assignment */
  for(direction=0;direction<3;direction++) {       
    assign_forces(1.0/(Mesh*Mesh*Mesh), F_K[direction], Teilchenzahl, Fmesh[direction],0);
    }

  return;
}
