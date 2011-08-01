
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>

#include "p3m.h"
#include "charge-assign.h"

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

void forward_fft(void) {
  fftw_execute(forward_plan);
}

void backward_fft(void) {
  int i;
  for(i=0;i<3;i++)
    fftw_execute(backward_plan[i]);
}



static FLOAT_TYPE sinc(FLOAT_TYPE d)
{
  /* 
     Berechnet die sinc-Funktion als sin(PI*x)/(PI*x).
     (Konvention fuer sinc wie in Hockney/Eastwood!)
  */

  static FLOAT_TYPE epsi = 1e-8;
  FLOAT_TYPE PId = PI*d;
  
  return (fabs(d)<=epsi) ? 1.0 : sin(PId)/PId;
}

void Init_interlaced_ik(int Teilchenzahl) {
  int l;
  Qmesh = (FLOAT_TYPE *) realloc(Qmesh, 2*Mesh*Mesh*Mesh*sizeof(FLOAT_TYPE));

  G_hat = (FLOAT_TYPE *) realloc(G_hat, Mesh*Mesh*Mesh*sizeof(G_hat));  

  for(l=0;l<2;l++) {
      ca_ind[l] = (int *) realloc(ca_ind[l], 3*Teilchenzahl*sizeof(int));
      cf[l] = (FLOAT_TYPE *) realloc(cf[l], Teilchenzahl * (ip+1) * (ip + 1) * (ip + 1) *sizeof(FLOAT_TYPE)); 
  }

  F_K[0] = Fx_K;
  F_K[1] = Fy_K;
  F_K[2] = Fz_K;

  forward_plan = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)Qmesh, (fftw_complex *)Qmesh, FFTW_FORWARD, FFTW_ESTIMATE);
  for(l=0;l<3;l++){
    Fmesh[l] = (FLOAT_TYPE *)realloc(Fmesh[l], 2*Mesh*Mesh*Mesh*sizeof(FLOAT_TYPE));  
    backward_plan[l] = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)Fmesh[l], (fftw_complex *)Fmesh[l], FFTW_BACKWARD, FFTW_ESTIMATE);
  }
  
}

void Aliasing_sums_interlaced_ik(int NX, int NY, int NZ, FLOAT_TYPE alpha,
				  FLOAT_TYPE *Zaehler, FLOAT_TYPE *Nenner1, FLOAT_TYPE *Nenner2)
{
  /*
    Berechnet die beiden Aliasing-Summen im Zaehler und Nenner des Ausdrucks fuer die
    optimale influence-function (siehe: Hockney/Eastwood: Formel 8-22 Seite 275 oben).
    
    NX,NY,NZ : Komponenten des n-Vektors, fuer den die Aliasing-Summen ausgefuehrt werden sollen.
    *ZaehlerX,*ZaehlerY,*ZaehlerZ : x- ,y- und z-Komponente der Aliasing-Summe im Zaehler.
    *Nenner : Aliasing-Summe im Nenner.
  */
  
  static int aliasmax = 1; /* Genauigkeit der Aliasing-Summe (2 ist wohl genug) */
  
  FLOAT_TYPE S,S1,S2,S3;
  FLOAT_TYPE fak1,fak2,zwi;
  int    MX,MY,MZ;
  FLOAT_TYPE NMX,NMY,NMZ;
  FLOAT_TYPE NM2;

  fak1 = 1.0/(FLOAT_TYPE)Mesh;
  fak2 = SQR(PI/(alpha*Len));

  Zaehler[0] = Zaehler[1] = Zaehler[2] = *Nenner1 = *Nenner2 = 0.0;

  for (MX = -aliasmax; MX <= aliasmax; MX++)
    {
      NMX = nshift[NX] + Mesh*MX;
      S   = SQR(sinc(fak1*NMX)); 
      S1  = pow(S,ip+1);
      for (MY = -aliasmax; MY <= aliasmax; MY++)
	{
	  NMY = nshift[NY] + Mesh*MY;
	  S   = SQR(sinc(fak1*NMY));
	  S2  = S1*pow(S,ip+1);
	  for (MZ = -aliasmax; MZ <= aliasmax; MZ++)
	    {
	      NMZ = nshift[NZ] + Mesh*MZ;
	      S   = SQR(sinc(fak1*NMZ));
	      S3  = S2*pow(S,ip+1);
	      NM2 = SQR(NMX) + SQR(NMY) + SQR(NMZ);

              *Nenner1 += S3;

	      zwi  = S3 * exp(-fak2*NM2)/NM2;
	      Zaehler[0] += NMX*zwi;
              Zaehler[1] += NMY*zwi;
              Zaehler[2] += NMZ*zwi;
	      
	       if (((MX+MY+MZ)%2)==0) {					//even term
	        *Nenner2 += S3;
	       } else {						//odd term: minus sign!
	         *Nenner2 -= S3;
	       }
	    }
	}
    }
}

 
void Influence_function_berechnen_ik_interlaced(FLOAT_TYPE alpha)
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
  FLOAT_TYPE Zaehler[3]={0.0,0.0,0.0},Nenner1=0.0, Nenner2=0.0;
  FLOAT_TYPE zwi;
  
  dMesh = (FLOAT_TYPE)Mesh;
  dMeshi= 1.0/dMesh;
  
  /* bei Zahlen >= Mesh/2 wird noch Mesh abgezogen! */
  for (NX=0; NX<Mesh; NX++)
    {
      for (NY=0; NY<Mesh; NY++)
	{
	  for (NZ=0; NZ<Mesh; NZ++)
	    {
	      if ((NX==0) && (NY==0) && (NZ==0))
		G_hat[r_ind(NX,NY,NZ)]=0.0;
              else if ((NX%(Mesh/2) == 0) && (NY%(Mesh/2) == 0) && (NZ%(Mesh/2) == 0))
                G_hat[r_ind(NX,NY,NZ)]=0.0;
	      else
		{
		  Aliasing_sums_interlaced_ik(NX,NY,NZ,alpha,Zaehler,&Nenner1, &Nenner2);
		  
		  Dnx = Dn[NX];  
		  Dny = Dn[NY];  
		  Dnz = Dn[NZ];
		  
		  zwi  = Dnx*Zaehler[0] + Dny*Zaehler[1] + Dnz*Zaehler[2];
		  zwi /= ( SQR(Dnx) + SQR(Dny) + SQR(Dnz) );
                  zwi /= 0.5*(SQR(Nenner1) + SQR(Nenner2));		  

		  G_hat[r_ind(NX,NY,NZ)] = Mesh*Mesh*Mesh*2.0 * zwi * Leni * Leni;
		}
	    }
	}
    }
}


void P3M_ik_interlaced(const FLOAT_TYPE alpha, const int Teilchenzahl)
{
  /*
    Berechnet den k-Raum Anteil der Ewald-Routine auf dem Gitter.
    Die Ableitung wird durch analytische Differentiation der charge assignment
    function erreicht, so wie es auch im EPBDLP-paper geschieht.
    alpha : Ewald-Parameter.
  */
  
  /* Zaehlvariablen: */
  int i, j, k, l, m, ii; 
  /* wahre n-Werte im Impulsraum: */
  int nx,ny,nz;
  /* Variablen fuer FFT: */
  int Nx,Ny,Nz,Lda, Ldb, sx, sy, sz, status;   
  /* Schnelles Modulo: */
  int MESHMASKE;
  /* Hilfsvariablen */
  FLOAT_TYPE dx,dy,dz,d2,H,Hi,dMesh,MI2,modadd2,modadd3;
  int modadd1;
  /* charge-assignment beschleunigen */
  FLOAT_TYPE T1,T2,T3,T4,T5;
  /* schnellerer Zugriff auf die Arrays Gx[i] etc.: */
  int Gxi,Gyi,Gzi;
  /* Argumente fuer das Array LadInt */
  int xarg,yarg,zarg;
  /* Gitterpositionen */
  int xpos,ypos,zpos;
  /* Soweit links vom Referenzpunkt gehts beim Ladungsver-
     teilen los (implementiert ist noch ein Summand Mesh!): */
  int mshift;
  FLOAT_TYPE dTeilchenzahli;
  
  //Pour le graphique des selfforces:
  int point,nb_points;
  FLOAT_TYPE posx=0,posy=0,posz=0;

  FILE *frho;

  Nx = Ny = Nz = Mesh;
  Lda = Ldb = Mesh; 
  sx = sy = sz = 1; 
  MESHMASKE = Mesh-1;
  dMesh = (FLOAT_TYPE)Mesh;
  H = Len/dMesh;
  Hi = 1.0/H;
  MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;
  mshift = Mesh-ip/2;
  dTeilchenzahli = 1.0/(FLOAT_TYPE)Teilchenzahl;

  /* Initialisieren von Q_re und Q_im */
  for(i=0;i<(2*Mesh*Mesh*Mesh); i++)
    Qmesh[i] = 0.0;

  /* chargeassignment */
  for (i=0; i<Teilchenzahl; i++)
    {
      for(ii=0;ii<2;ii++) {
        FLOAT_TYPE rpos[3] = {xS[i] + 0.5*ii*H, yS[i] + 0.5*ii*H, zS[i] + 0.5*ii*H};
	assign_charge(i, Q[i], rpos, Qmesh, ii);
      }
    }
  
  /* Durchfuehren der Fourier-Hin-Transformationen: */
  forward_fft();

  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<Mesh; k++)
	{
          int c_index = c_ind(i,j,k);
          FLOAT_TYPE dop;
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
   
  /* Kraftkomponenten: */
  for(i=0;i<3;i++) { 
    assign_forces( 1.0/(Mesh*Mesh*Mesh), F_K[i], Teilchenzahl, Fmesh[i], 0);
    assign_forces( 1.0/(Mesh*Mesh*Mesh), F_K[i], Teilchenzahl, Fmesh[i], 1);
  }
  return;
}
