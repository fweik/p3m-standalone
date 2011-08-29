
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>

#include "p3m.h"
#include "common.h"
#include "p3m-common.h"
#include "charge-assign.h"
#include "p3m-ik-i.h"

// declaration of the method

const method_t method_p3m_ik_i = { MEHOTD_P3M_ik_i, "P3M with ik differentiation, interlaced.",
				   METHOD_FLAG_ik | METHOD_FLAG_Qmesh | METHOF_FLAG_G_hat | METHOD_FLAG_nshift | METHOD_FLAG_ca | METHOD_FLAG_interlaced,
				   &Init_ik_i, &Influence_function_ik_i, &P3M_ik_i, NULL };

fftw_plan forward_plan;
fftw_plan backward_plan[3];

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

void Init_interlaced_ik( system_t *s, parameters_t *p ) {
  int l;
  int Mesh = p->mesh;

  data_t *d = Init_data( &method_p3m_ik_i, s, p );

  forward_plan = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)d->Qmesh, (fftw_complex *)d->Qmesh, FFTW_FORWARD, FFTW_ESTIMATE);
  for(l=0;l<3;l++){
    backward_plan[l] = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)d->Fmesh.fields[l], (fftw_complex *)d->Fmesh.fields[l], FFTW_BACKWARD, FFTW_ESTIMATE);
  }
  
}

void Aliasing_sums_interlaced_ik( system_t *s, parameters_t *p, data_t *d, int NX, int NY, int NZ,
				  FLOAT_TYPE *Zaehler, FLOAT_TYPE *Nenner1, FLOAT_TYPE *Nenner2)
{
  FLOAT_TYPE S,S1,S2,S3;
  FLOAT_TYPE fak1,fak2,zwi;
  int    MX,MY,MZ;
  FLOAT_TYPE NMX,NMY,NMZ;
  FLOAT_TYPE NM2;
  
  FLOAT_TYPE expo, TE;
  int Mesh = d->mesh;
  FLOAT_TYPE Leni = 1.0/s->length;

  fak1 = 1.0/(FLOAT_TYPE)Mesh;
  fak2 = SQR(PI/(alpha*Len));
  
  Zaehler[0] = Zaehler[1] = Zaehler[2] = *Nenner1 = *Nenner2 = 0.0;

  for (MX = -P3M_BRILLOUIN; MX <= P3M_BRILLOUIN; MX++)
    {
      NMX = nshift[NX] + Mesh*MX;
      S   = SQR(sinc(fak1*NMX)); 
      S1  = pow(S,ip+1);
      for (MY = -P3M_BRILLOUIN; MY <= P3M_BRILLOUIN; MY++)
	{
	  NMY = nshift[NY] + Mesh*MY;
	  S   = SQR(sinc(fak1*NMY));
	  S2  = S1*pow(S,ip+1);
	  for (MZ = -P3M_BRILLOUIN; MZ <= P3M_BRILLOUIN; MZ++)
	    {
	      NMZ = nshift[NZ] + Mesh*MZ;
	      S   = SQR(sinc(fak1*NMZ));
	      S3  = S2*pow(S,ip+1);
	      NM2 = SQR(NMX) + SQR(NMY) + SQR(NMZ);

              *Nenner1 += S3;

	      expo = fak2*NM2;
              TE = ( expo < 30.0 ) ? exp ( -expo ) : 0.0;
	      zwi  = S3 * TE/NM2;
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

 
void Influence_function_berechnen_ik_interlaced( system_t *s, parameters_t *p, data_t *d )
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

  int ind = 0;
  int Mesh = p->mesh;
  FLOAT_TYPE Leni = 1.0/s->length;  

  dMesh = (FLOAT_TYPE)Mesh;
  dMeshi= 1.0/dMesh;
  
  /* bei Zahlen >= Mesh/2 wird noch Mesh abgezogen! */
  for (NX=0; NX<Mesh; NX++)
    {
      for (NY=0; NY<Mesh; NY++)
	{
	  for (NZ=0; NZ<Mesh; NZ++)
	    {
	      ind = r_ind( NX, NZ, NZ );

	      if ((NX==0) && (NY==0) && (NZ==0))
		d->G_hat[ind]=0.0;
              else if ((NX%(Mesh/2) == 0) && (NY%(Mesh/2) == 0) && (NZ%(Mesh/2) == 0))
                d->G_hat[ind]=0.0;
	      else
		{
		  Aliasing_sums_interlaced_ik( s, p, d, NX, NY, NZ, Zaehler, &Nenner1,  &Nenner2);
		  
		  Dnx = Dn[NX];  
		  Dny = Dn[NY];  
		  Dnz = Dn[NZ];
		  
		  zwi  = Dnx*Zaehler[0] + Dny*Zaehler[1] + Dnz*Zaehler[2];
		  zwi /= ( SQR(Dnx) + SQR(Dny) + SQR(Dnz) );
                  zwi /= 0.5*(SQR(Nenner1) + SQR(Nenner2));		  
		  
		  G_hat[ind] = Mesh*Mesh*Mesh*2.0 * zwi * Leni * Leni;
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
