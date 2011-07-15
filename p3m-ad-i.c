#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include <fftw3.h>

#include "p3m.h"
#include "p3m-ad-i.h"

//VB: define flag INTERLACED to activate interlacing:
//VB: define flag SUBSTRACT_SF to substract self-forces:

#define SUBSTRACT_SF

/* Pi, weil man's so oft braucht: */
#define PI 3.14159265358979323846264

fftw_plan forward_plan;
fftw_plan backward_plan;

FLOAT_TYPE Potenz(FLOAT_TYPE x, int ip);
FLOAT_TYPE sinc(FLOAT_TYPE d);

FLOAT_TYPE *Q_re, *Q_im;

//VB: calcul des self-forces
#ifdef SUBSTRACT_SF
FLOAT_TYPE *SelfForce_x,*SelfForce_y,*SelfForce_z;
#endif

//VB: for interlacing:
FLOAT_TYPE *F2x_K, *F2y_K, *F2z_K;

#define c_idx(A,B,C) (2*Mesh*Mesh*(A) + 2*Mesh * (B) + 2*(C))
#define r_ind(A,B,C) ((A)*Mesh*Mesh + (B)*Mesh + (C))

void forward_fft(void) {
  int i,j,k;
  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<Mesh; k++) {
        Qmesh[c_idx(i,j,k)] = Q_re[r_ind(i,j,k)];
        Qmesh[c_idx(i,j,k)+1] = 0.0;
      }
  fftw_execute(forward_plan);

  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<Mesh; k++) {
        Q_re[r_ind(i,j,k)] = Qmesh[c_idx(i,j,k)];
        Q_im[r_ind(i,j,k)] = Qmesh[c_idx(i,j,k)+1];
      }
}

void backward_fft(void) {
  int i,j,k;
  
  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<Mesh; k++) {
        Qmesh[c_idx(i,j,k)] = Q_re[r_ind(i,j,k)];
        Qmesh[c_idx(i,j,k)+1] = Q_im[r_ind(i,j,k)];
      }

  fftw_execute(backward_plan);

  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<Mesh; k++) {
        Q_re[r_ind(i,j,k)] = Qmesh[c_idx(i,j,k)];
        Q_im[r_ind(i,j,k)] = Qmesh[c_idx(i,j,k)+1];
        if(Q_im[r_ind(i,j,k)] > 10e-5)
          puts("Warning, non-real potential.");
      }
}


//VB

void Aliasing_sums_AD_interlacing(int NX, int NY, int NZ, FLOAT_TYPE alpha,
		   FLOAT_TYPE *Zaehler, FLOAT_TYPE *Nenner1, FLOAT_TYPE *Nenner2, FLOAT_TYPE *Nenner3, FLOAT_TYPE *Nenner4)
{
  /*
    Berechnet die beiden Aliasing-Summen im Zaehler und Nenner des Ausdrucks fuer die
    optimale influence-function (siehe: Hockney/Eastwood: Formel 8-22 Seite 275 oben).
    
    NX,NY,NZ : Komponenten des n-Vektors, fuer den die Aliasing-Summen ausgefuehrt werden sollen.
    *ZaehlerX,*ZaehlerY,*ZaehlerZ : x- ,y- und z-Komponente der Aliasing-Summe im Zaehler.
    *Nenner : Aliasing-Summe im Nenner.
  */
  
  static int aliasmax = 2; /* Genauigkeit der Aliasing-Summe (2 ist wohl genug) */
  
  FLOAT_TYPE S,S1,S2,S3;
  FLOAT_TYPE fak1,fak2,zwi;
  int    MX,MY,MZ;
  FLOAT_TYPE NMX,NMY,NMZ;
  FLOAT_TYPE NM2;

  fak1 = 1.0/(FLOAT_TYPE)Mesh;
  fak2 = SQR(PI/(alpha*Len));

  *Zaehler = *Nenner1 = *Nenner2 = *Nenner3 = *Nenner4 = 0.0;

  for (MX = -aliasmax; MX <= aliasmax; MX++)
    {
      NMX = nshift[NX] + Mesh*MX;
      S   = SQR(sinc(fak1*NMX)); 
      S1  = Potenz(S,ip);
      for (MY = -aliasmax; MY <= aliasmax; MY++)
	{
	  NMY = nshift[NY] + Mesh*MY;
	  S   = SQR(sinc(fak1*NMY));
	  S2  = S1*Potenz(S,ip);
	  for (MZ = -aliasmax; MZ <= aliasmax; MZ++)
	    {

	      NMZ = nshift[NZ] + Mesh*MZ;
	      S   = SQR(sinc(fak1*NMZ));
	      S3  = S2*Potenz(S,ip);
	      NM2 = SQR(NMX) + SQR(NMY) + SQR(NMZ);
	      
	      *Nenner1 += S3;
	      *Nenner2 += S3 * NM2;
	      
	      zwi  = S3 * exp(-fak2*NM2);
	      *Zaehler += zwi;
	      
	 if (((MX+MY+MZ)%2)==0) {					//even term
	   *Nenner3 += S3;
	   *Nenner4 += S3*NM2;
	 } else {						//odd term: minus sign!
	   *Nenner3 -= S3;
	   *Nenner4 -= S3*NM2;
	 }


	    }
	}
    }
}

/*------------------------------------------------------------*/

void Influence_function_ad_interlaced(FLOAT_TYPE alpha)
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
  FLOAT_TYPE fak1,fak2,dMesh,dMeshi;
  FLOAT_TYPE Zaehler,Nenner1,Nenner2,Nenner3,Nenner4;
  FLOAT_TYPE zwi;
  
  FLOAT_TYPE qua,qua_;
  
  dMesh = (FLOAT_TYPE)Mesh;
  dMeshi= 1.0/dMesh;
  
  fak1 = 2.0/(SQR(Len) * 2.0*PI);
  fak2 = SQR(PI/(alpha*Len));
  
  /* bei Zahlen >= Mesh/2 wird noch Mesh abgezogen! */
  for (NX=0; NX<Mesh; NX++)
    {
      fprintf(stderr,"."); fflush(stderr);
      for (NY=0; NY<Mesh; NY++)
	{
	  for (NZ=0; NZ<Mesh; NZ++)
	    {
	      if ((NX==0) && (NY==0) && (NZ==0))
		G_hat[r_ind(NX,NY,NZ)]=0.0;
	      else
		{
		  Aliasing_sums_AD_interlacing(NX,NY,NZ,alpha,&Zaehler,&Nenner1,&Nenner2,&Nenner3,&Nenner4);
		  Dnx = Dn[NX];  
		  Dny = Dn[NY];  
		  Dnz = Dn[NZ];
		  
		  zwi  = Zaehler;
		  zwi /= ( 0.5*(Nenner1*Nenner2 + Nenner3*Nenner4) );
		  
		  G_hat[r_ind(NX,NY,NZ)] = zwi*fak1;
		}
	    }
	}
    }
  fprintf(stderr,"\n");
}


void P3M_ad_interlaced(FLOAT_TYPE alpha, int Teilchenzahl)
{
  /*
    Berechnet den k-Raum Anteil der Ewald-Routine auf dem Gitter.
    Die Ableitung wird durch analytische Differentiation der charge assignment
    function erreicht, so wie es auch im EPBDLP-paper geschieht.
    alpha : Ewald-Parameter.
  */
  
  /* Zaehlvariablen: */
  int i, j, k, l, m; 
  /* wahre n-Werte im Impulsraum: */
  int nx,ny,nz;
  /* Fouriertransformierte Differenzenoperatoren */
  FLOAT_TYPE Dnx,Dny,Dnz;
  /* Variablen fuer FFT: */
  int sx, sy, sz, status;   
  /* Schnelles Modulo: */
  int MESHMASKE;
  /* Hilfsvariablen */
  FLOAT_TYPE d1,d2,H,Hi,dMesh,MI2,modadd1,modadd2,modadd3;
  /* charge-assignment beschleunigen */
  FLOAT_TYPE T1,T2,T3,T4,T5;
  FLOAT_TYPE T1_x,T2_x,T3_x;
  FLOAT_TYPE T1_y,T2_y,T3_y;
  FLOAT_TYPE T1_z,T2_z,T3_z;
  /* schnellerer Zugriff auf die Arrays Gx[i] etc.: */
  int Gxi,Gyi,Gzi;
  /* Argumente fuer das Array LadInt */
  int xarg,yarg,zarg;
  /* Gitterpositionen */
  int xpos,ypos,zpos;
  /* Soweit links vom Referenzpunkt gehts beim Ladungsver-
     teilen los (implementiert ist noch ein Summand Mesh!): */
  int mshift;
  FLOAT_TYPE Energiefaktor,Kraftfaktor,dTeilchenzahli;
  
  //VB: Variables pour les self-forces
  int ii;
  
  //Pour le graphique des selfforces:
  int point,nb_points;
  FLOAT_TYPE posx=0,posy=0,posz=0;
  FILE* fout;

  sx = sy = sz = 1; 
  MESHMASKE = Mesh-1;
  dMesh = (FLOAT_TYPE)Mesh;
  H = Len/dMesh;
  Hi = 1.0/H;
  MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;
  mshift = Mesh-ip/2;
  dTeilchenzahli = 1.0/(FLOAT_TYPE)Teilchenzahl;
  Energiefaktor  = 0.5*Len/(dMesh*dMesh*dMesh);
  
    /* Vorbereitung der Fallunterscheidung nach geradem/ungeradem ip: */
  switch (ip)
    {
    case 0 : case 2 : case 4 : case 6 :
      { modadd1=0.0; modadd2=0.5; modadd3= 0.5;} break;
    case 1 :case 3 : case 5 :
      { modadd1=0.5; modadd2=0.0; modadd3=-0.5;} break;
      /* Beachte: modadd1+modadd3=modadd2! */
    default : 
      {
	fprintf(stderr,"Wert von ip (ip=%d) ist nicht erlaubt!",ip);
	fprintf(stderr,"Programm abgebrochen!");
	exit(1);
      } break;
    }
    
//============================================
//
//      Premier calcul des forces
//
//============================================

  

//VB:
  //--------------------------
  // Graphique des self-forces
  //--------------------------
#ifdef SUBSTRACT_SF
  printf("Calcul self-forces\n");
  //--------------------------
  // Calcul des self-forces
  //--------------------------
  for (ii=0; ii<Teilchenzahl; ii++) {
//Print force:
//printf("AVANT: self-force[%i] = %1.9lf %1.9lf %1.9lf\n",ii,SelfForce_x[ii],SelfForce_y[ii],SelfForce_z[ii]);

     m = 0;
       // Initialisieren von Q_re und Q_im 
     for (i=0; i<Mesh; i++)
       for (j=0; j<Mesh; j++)
	 for (k=0; k<Mesh; k++)
	   Q_re[r_ind(i,j,k)] = Q_im[r_ind(i,j,k)] = 0.0;



      d1 = xS[ii]*Hi + modadd1; 
      Gx[ii] = Gxi = (int)(d1 + modadd3) + mshift;
      xarg = (int)( (d1 - round(d1) + 0.5)*MI2 );
      
      d1 = yS[ii]*Hi + modadd1; 
      Gy[ii] = Gyi = (int)(d1 + modadd3) + mshift;
      yarg = (int)( (d1 - round(d1) + 0.5)*MI2 );
      
      d1 = zS[ii]*Hi + modadd1; 
      Gz[ii] = Gzi = (int)(d1 + modadd3) + mshift;
      zarg = (int)( (d1 - round(d1) + 0.5)*MI2 );
      
      for (j=0; j<=ip; j++)
	{
	  xpos = (Gxi+j)&MESHMASKE;
	  T1   = Q[ii]*LadInt[j][xarg];
	  T1_x = dMesh * Q[ii] * LadInt_[j][xarg];
	  T1_y = dMesh * T1; 
	  T1_z = dMesh * T1; 
	  for (k=0; k<=ip; k++)
	    {
	      ypos = (Gyi+k)&MESHMASKE;
	      T4   = LadInt[k][yarg];
	      T2   = T1   * T4;
	      T2_x = T1_x * T4;
	      T2_y = T1_y * LadInt_[k][yarg];
	      T2_z = T1_z * T4;
	      for (l=0; l<=ip; l++)
		{
		  zpos = (Gzi+l)&MESHMASKE;
		  T5   = LadInt[l][zarg];
		  T3   = T2   * T5;
		  T3_x = T2_x * T5;
		  T3_y = T2_y * T5;
		  T3_z = T2_z * LadInt_[l][zarg];

		  Q_re[r_ind(xpos,ypos,zpos)] += T3;
		  dQdx[m] = T3_x;
		  dQdy[m] = T3_y;
		  dQdz[m] = T3_z;

		  m++;
		}
	    }
	}
      
        // Durchfuehren der Fourier-Hin-Transformationen: 
   
      forward_fft();
       
  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<Mesh; k++)
	{
	  T1 = G_hat[r_ind(i,j,k)];
	  Q_re[r_ind(i,j,k)] *= T1;
	  Q_im[r_ind(i,j,k)] *= T1;
	}

  backward_fft();

  
  // Durchfuehren der Fourier-Rueck-Transformation: 
   
  // Kraftkomponenten:
  m = 0;
	 for (j=0; j<=ip; j++) 
	 {
	  xpos = (Gx[ii]+j)&MESHMASKE;
	  for (k=0; k<=ip; k++) 
	    {
	      ypos = (Gy[ii]+k)&MESHMASKE;
	      for (l=0; l<=ip; l++) 
		{
		  zpos = (Gz[ii]+l)&MESHMASKE;

		  T1 = Q_re[r_ind(xpos,ypos,zpos)];

		  SelfForce_x[ii] -= T1*dQdx[m];
		  SelfForce_y[ii] -= T1*dQdy[m];
		  SelfForce_z[ii] -= T1*dQdz[m];

		  m++;
		}
	    }
	 }
     //Print force:
//     printf("self-force[%i] = %1.9lf %1.9lf %1.9lf\n",ii,SelfForce_x[ii],SelfForce_y[ii],SelfForce_z[ii]);
  }
  //Fin des self-forces
  //------------------------------------------------------------------------
#endif

  puts("Assigning charges...");

  /* Verteilung der Ladungen auf's Gitter: */
  
  /* Initialisieren von Q_re und Q_im */
  for (i=0; i<(Mesh*Mesh*Mesh); i++)
	Q_re[i] = Q_im[i] = 0.0;


  /* Die Arrays Q, dQdx, dQdy und dQdz berechnen */
  m = 0;   /* <-- Zaehlt alle fraktionellen Ladungsbeitraege linear durch! */
  for (i=0; i<Teilchenzahl; i++)
    {
      d1 = xS[i]*Hi + modadd1; 
      Gx[i] = Gxi = (int)(d1 + modadd3) + mshift;
      xarg = (int)( (d1 - round(d1) + 0.5)*MI2 );
      
      d1 = yS[i]*Hi + modadd1; 
      Gy[i] = Gyi = (int)(d1 + modadd3) + mshift;
      yarg = (int)( (d1 - round(d1) + 0.5)*MI2 );
      
      d1 = zS[i]*Hi + modadd1; 
      Gz[i] = Gzi = (int)(d1 + modadd3) + mshift;
      zarg = (int)( (d1 - round(d1) + 0.5)*MI2 );
      
      for (j=0; j<=ip; j++)
	{
	  xpos = (Gxi+j)&MESHMASKE;
	  T1   = Q[i]*LadInt[j][xarg];
	  T1_x = dMesh * Q[i] * LadInt_[j][xarg];
	  T1_y = dMesh * T1; 
	  T1_z = dMesh * T1; 
	  for (k=0; k<=ip; k++)
	    {
	      ypos = (Gyi+k)&MESHMASKE;
	      T4   = LadInt[k][yarg];
	      T2   = T1   * T4;
	      T2_x = T1_x * T4;
	      T2_y = T1_y * LadInt_[k][yarg];
	      T2_z = T1_z * T4;
	      for (l=0; l<=ip; l++)
		{
		  zpos = (Gzi+l)&MESHMASKE;
		  T5   = LadInt[l][zarg];
		  T3   = T2   * T5;
		  T3_x = T2_x * T5;
		  T3_y = T2_y * T5;
		  T3_z = T2_z * LadInt_[l][zarg];

		  Q_re[r_ind(xpos,ypos,zpos)] += T3;
		  dQdx[m] = T3_x;
		  dQdy[m] = T3_y;
		  dQdz[m] = T3_z;

		  m++;
		}
	    }
	}
    }
  puts("Finished charge assignment.");
  /* Durchfuehren der Fourier-Hin-Transformationen: */

  forward_fft();
   
  puts("forward fft done.");

  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<Mesh; k++)
	{
	  T1 = G_hat[r_ind(i,j,k)];
	  E_Coulomb_Impuls_Summe += T1*(SQR(Q_re[r_ind(i,j,k)]) + SQR(Q_im[r_ind(i,j,k)]));
	  Q_re[r_ind(i,j,k)] *= T1;
	  Q_im[r_ind(i,j,k)] *= T1;
	}
  E_Coulomb_Impuls_Summe *= Energiefaktor;
  
  /* Durchfuehren der Fourier-Rueck-Transformation: */

  puts("kspace caclculation done");

  backward_fft();
   
  puts("backward fft done."); 

  /* Kraftkomponenten: */
  m = 0;
  for (i=0; i<Teilchenzahl; i++)
    {
      for (j=0; j<=ip; j++) 
	{
	  xpos = (Gx[i]+j)&MESHMASKE;
	  for (k=0; k<=ip; k++) 
	    {
	      ypos = (Gy[i]+k)&MESHMASKE;
	      for (l=0; l<=ip; l++) 
		{
		  zpos = (Gz[i]+l)&MESHMASKE;
		  
		  T1 = Q_re[r_ind(xpos,ypos,zpos)];
		  
		  Fx_K[i] -= T1*dQdx[m];
		  Fy_K[i] -= T1*dQdy[m];
		  Fz_K[i] -= T1*dQdz[m];
		  
		  m++;
		}
	    }
	}
    }
    
  puts("Force assignment done.");
//  for (i=0; i<Teilchenzahl; i++) {
//   printf("k-force[%i] = %1.9lf %1.9lf %1.9lf\n",i,Fx_K[i],Fy_K[i],Fz_K[i]);
//  }

  
#ifdef SUBSTRACT_SF
printf("Substract SF\n");
   // VB: SUBSTRACT SELF-FORCES:
 //printf("After substraction of self-force:\n");
  for (i=0; i<Teilchenzahl; i++) {
     Fx_K[i] -= SelfForce_x[i];
     Fy_K[i] -= SelfForce_y[i];
     Fz_K[i] -= SelfForce_z[i];
//     printf("k-force[%i] = %f %f %f\n",i,Fx_K[i],Fy_K[i],Fz_K[i]);
  }
#endif

  /* Kraftdrift beseitigen */
  T1 = T2 = T3 = 0.0;
  for (i=0; i<Teilchenzahl; i++) { T1 += Fx_K[i]; T2 += Fy_K[i]; T3 += Fz_K[i]; }
  

  //Force resultante:
//  printf("La resultante des forces est: %e %e %e\n",T1, T2, T3);
  
  T1 *= dTeilchenzahli; T2 *= dTeilchenzahli; T3 *= dTeilchenzahli;
  for (i=0; i<Teilchenzahl; i++) {
     Fx_K[i] -= T1; Fy_K[i] -= T2; Fz_K[i] -= T3;
  }

//========================================================================
//
//      Deuxieme calcul des forces (avec shift de h/2 selon la diagonale)
//
//========================================================================

#ifdef SUBSTRACT_SF
//VB: re-initialiser les self-forces:
  for (i=0; i<Teilchenzahl; i++)   {
      SelfForce_x[i] = 0.0;
      SelfForce_y[i] = 0.0;
      SelfForce_z[i] = 0.0;
   }


  //--------------------------
  // Calcul des self-forces
  //--------------------------
  for (ii=0; ii<Teilchenzahl; ii++) {
//Print force:
//printf("AVANT: self-force[%i] = %1.9lf %1.9lf %1.9lf\n",ii,SelfForce_x[ii],SelfForce_y[ii],SelfForce_z[ii]);

     m = 0;
       // Initialisieren von Q_re und Q_im 
     for (i=0; i<Mesh; i++)
       for (j=0; j<Mesh; j++)
	 for (k=0; k<Mesh; k++)
	   Q_re[r_ind(i,j,k)] = Q_im[r_ind(i,j,k)] = 0.0;


      d1 = (xS[ii]+H/2.0)*Hi + modadd1; 		//Shift de h/2 selon la diagonale!
      Gx[ii] = Gxi = (int)(d1 + modadd3) + mshift;
      xarg = (int)( (d1 - round(d1) + 0.5)*MI2 );
      
      d1 = (yS[ii]+H/2.0)*Hi + modadd1; 
      Gy[ii] = Gyi = (int)(d1 + modadd3) + mshift;
      yarg = (int)( (d1 - round(d1) + 0.5)*MI2 );
      
      d1 = (zS[ii]+H/2.0)*Hi + modadd1; 
      Gz[ii] = Gzi = (int)(d1 + modadd3) + mshift;
      zarg = (int)( (d1 - round(d1) + 0.5)*MI2 );
      
      for (j=0; j<=ip; j++)
	{
	  xpos = (Gxi+j)&MESHMASKE;
	  T1   = Q[ii]*LadInt[j][xarg];
	  T1_x = dMesh * Q[ii] * LadInt_[j][xarg];
	  T1_y = dMesh * T1; 
	  T1_z = dMesh * T1; 
	  for (k=0; k<=ip; k++)
	    {
	      ypos = (Gyi+k)&MESHMASKE;
	      T4   = LadInt[k][yarg];
	      T2   = T1   * T4;
	      T2_x = T1_x * T4;
	      T2_y = T1_y * LadInt_[k][yarg];
	      T2_z = T1_z * T4;
	      for (l=0; l<=ip; l++)
		{
		  zpos = (Gzi+l)&MESHMASKE;
		  T5   = LadInt[l][zarg];
		  T3   = T2   * T5;
		  T3_x = T2_x * T5;
		  T3_y = T2_y * T5;
		  T3_z = T2_z * LadInt_[l][zarg];

		  Q_re[r_ind(xpos,ypos,zpos)] += T3;
		  dQdx[m] = T3_x;
		  dQdy[m] = T3_y;
		  dQdz[m] = T3_z;

		  m++;
		}
	    }
	}
      
        // Durchfuehren der Fourier-Hin-Transformationen: 
      forward_fft();

  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<Mesh; k++)
	{
	  T1 = G_hat[r_ind(i,j,k)];
	  Q_re[r_ind(i,j,k)] *= T1;
	  Q_im[r_ind(i,j,k)] *= T1;
	}
  
  // Durchfuehren der Fourier-Rueck-Transformation: 
  backward_fft();
   
  // Kraftkomponenten:
  m = 0;
	 for (j=0; j<=ip; j++) 
	 {
	  xpos = (Gx[ii]+j)&MESHMASKE;
	  for (k=0; k<=ip; k++) 
	    {
	      ypos = (Gy[ii]+k)&MESHMASKE;
	      for (l=0; l<=ip; l++) 
		{
		  zpos = (Gz[ii]+l)&MESHMASKE;

		  T1 = Q_re[r_ind(xpos,ypos,zpos)];

		  SelfForce_x[ii] -= T1*dQdx[m];
		  SelfForce_y[ii] -= T1*dQdy[m];
		  SelfForce_z[ii] -= T1*dQdz[m];

		  m++;
		}
	    }
	 }
     //Print force:
//     printf("self-force[%i] = %1.9lf %1.9lf %1.9lf\n",ii,SelfForce_x[ii],SelfForce_y[ii],SelfForce_z[ii]);
  }
  //Fin des self-forces
  //------------------------------------------------------------------------
#endif

  /* Verteilung der Ladungen auf's Gitter: */
  
  /* Initialisieren von Q_re und Q_im */
  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<Mesh; k++)
	Q_re[r_ind(i,j,k)] = Q_im[r_ind(i,j,k)] = 0.0;


  /* Die Arrays Q, dQdx, dQdy und dQdz berechnen */
  m = 0;   /* <-- Zaehlt alle fraktionellen Ladungsbeitraege linear durch! */
  for (i=0; i<Teilchenzahl; i++)
    {
      d1 = (xS[i]+H/2.0)*Hi + modadd1; 
      Gx[i] = Gxi = (int)(d1 + modadd3) + mshift;
      xarg = (int)( (d1 - round(d1) + 0.5)*MI2 );
      
      d1 = (yS[i]+H/2.0)*Hi + modadd1; 
      Gy[i] = Gyi = (int)(d1 + modadd3) + mshift;
      yarg = (int)( (d1 - round(d1) + 0.5)*MI2 );
      
      d1 = (zS[i]+H/2.0)*Hi + modadd1; 
      Gz[i] = Gzi = (int)(d1 + modadd3) + mshift;
      zarg = (int)( (d1 - round(d1) + 0.5)*MI2 );
      
      for (j=0; j<=ip; j++)
	{
	  xpos = (Gxi+j)&MESHMASKE;
	  T1   = Q[i]*LadInt[j][xarg];
	  T1_x = dMesh * Q[i] * LadInt_[j][xarg];
	  T1_y = dMesh * T1; 
	  T1_z = dMesh * T1; 
	  for (k=0; k<=ip; k++)
	    {
	      ypos = (Gyi+k)&MESHMASKE;
	      T4   = LadInt[k][yarg];
	      T2   = T1   * T4;
	      T2_x = T1_x * T4;
	      T2_y = T1_y * LadInt_[k][yarg];
	      T2_z = T1_z * T4;
	      for (l=0; l<=ip; l++)
		{
		  zpos = (Gzi+l)&MESHMASKE;
		  T5   = LadInt[l][zarg];
		  T3   = T2   * T5;
		  T3_x = T2_x * T5;
		  T3_y = T2_y * T5;
		  T3_z = T2_z * LadInt_[l][zarg];

		  Q_re[r_ind(xpos,ypos,zpos)] += T3;
		  dQdx[m] = T3_x;
		  dQdy[m] = T3_y;
		  dQdz[m] = T3_z;

		  m++;
		}
	    }
	}
    }

  /* Durchfuehren der Fourier-Hin-Transformationen: */
  forward_fft();

  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<Mesh; k++)
	{
	  T1 = G_hat[r_ind(i,j,k)];
	  E_Coulomb_Impuls_Summe += T1*(SQR(Q_re[r_ind(i,j,k)]) + SQR(Q_im[r_ind(i,j,k)]));
	  Q_re[r_ind(i,j,k)] *= T1;
	  Q_im[r_ind(i,j,k)] *= T1;
	}
  E_Coulomb_Impuls_Summe *= Energiefaktor;
  
  /* Durchfuehren der Fourier-Rueck-Transformation: */
  backward_fft();
   
  /* Kraftkomponenten: */
  m = 0;
  for (i=0; i<Teilchenzahl; i++)
    {
      for (j=0; j<=ip; j++) 
	{
	  xpos = (Gx[i]+j)&MESHMASKE;
	  for (k=0; k<=ip; k++) 
	    {
	      ypos = (Gy[i]+k)&MESHMASKE;
	      for (l=0; l<=ip; l++) 
		{
		  zpos = (Gz[i]+l)&MESHMASKE;
		  
		  T1 = Q_re[r_ind(xpos,ypos,zpos)];
		  
		  F2x_K[i] -= T1*dQdx[m];
		  F2y_K[i] -= T1*dQdy[m];
		  F2z_K[i] -= T1*dQdz[m];
		  
		  m++;
		}
	    }
	}
    }
    
//  for (i=0; i<Teilchenzahl; i++) {
//   printf("k-force[%i] = %1.9lf %1.9lf %1.9lf\n",i,F2x_K[i],F2y_K[i],F2z_K[i]);
//  }


#ifdef SUBSTRACT_SF
   // VB: SUBSTRACT SELF-FORCES:
 //printf("After substraction of self-force:\n");
  for (i=0; i<Teilchenzahl; i++) {
     F2x_K[i] -= SelfForce_x[i];
     F2y_K[i] -= SelfForce_y[i];
     F2z_K[i] -= SelfForce_z[i];
//     printf("k-force[%i] = %f %f %f\n",i,F2x_K[i],F2y_K[i],F2z_K[i]);
  }
#endif

  /* Kraftdrift beseitigen */
  T1 = T2 = T3 = 0.0;
  for (i=0; i<Teilchenzahl; i++) { T1 += F2x_K[i]; T2 += F2y_K[i]; T3 += F2z_K[i]; }
  

  //Force resultante:
//  printf("La resultante des forces est: %e %e %e\n",T1, T2, T3);
  
  T1 *= dTeilchenzahli; T2 *= dTeilchenzahli; T3 *= dTeilchenzahli;
  for (i=0; i<Teilchenzahl; i++) {
     F2x_K[i] -= T1; F2y_K[i] -= T2; F2z_K[i] -= T3;
  }
  
//================= Calcul de la moyenne des deux forces =================

  for (i=0; i<Teilchenzahl; i++) {
     Fx_K[i] = (Fx_K[i]+F2x_K[i])/2.0;
     Fy_K[i] = (Fy_K[i]+F2y_K[i])/2.0;
     Fz_K[i] = (Fz_K[i]+F2z_K[i])/2.0;
//     printf("k-force[%i] = %f %f %f\n",i,Fx_K[i],Fy_K[i],Fz_K[i]);
  }
  
  return;
  
  //Not called;

}

/*------------------------------------------------------------*/

FLOAT_TYPE Potenz(FLOAT_TYPE x, int ip)
{
  /* Liefert x^(ip+1) zurueck */
  
  return pow(x, ip+1.0);
}

FLOAT_TYPE sinc(FLOAT_TYPE d)
{
  /* 
     Berechnet die sinc-Funktion als sin(PI*x)/(PI*x).
     (Konvention fuer sinc wie in Hockney/Eastwood!)
  */

  static FLOAT_TYPE epsi = 1e-8;
  FLOAT_TYPE PId = PI*d;
  
  return (fabs(d)<=epsi) ? 1.0 : sin(PId)/PId;
}

void Init_ad_interlaced(int Teilchenzahl) {
  int l = (ip+1)*(ip+1)*(ip+1);
  Qmesh = (FLOAT_TYPE *) realloc(Qmesh, 2*Mesh*Mesh*Mesh*sizeof(FLOAT_TYPE));
  Q_re = (FLOAT_TYPE *) realloc(Q_re, Mesh*Mesh*Mesh*sizeof(FLOAT_TYPE));
  Q_im = (FLOAT_TYPE *) realloc(Q_im, Mesh*Mesh*Mesh*sizeof(FLOAT_TYPE));
  Gx = (int *)realloc(Gx, Teilchenzahl*sizeof(int));
  Gy = (int *)realloc(Gy, Teilchenzahl*sizeof(int));
  Gz = (int *)realloc(Gz, Teilchenzahl*sizeof(int));

  dQdx = (FLOAT_TYPE *)realloc(dQdx, l*Teilchenzahl*sizeof(FLOAT_TYPE));
  dQdy = (FLOAT_TYPE *)realloc(dQdy, l*Teilchenzahl*sizeof(FLOAT_TYPE));
  dQdz = (FLOAT_TYPE *)realloc(dQdz, l*Teilchenzahl*sizeof(FLOAT_TYPE));

  F2x_K = (FLOAT_TYPE *) realloc(F2x_K, Teilchenzahl*sizeof(FLOAT_TYPE));
  F2y_K = (FLOAT_TYPE *) realloc(F2y_K, Teilchenzahl*sizeof(FLOAT_TYPE));
  F2z_K = (FLOAT_TYPE *) realloc(F2z_K, Teilchenzahl*sizeof(FLOAT_TYPE));

  G_hat = (FLOAT_TYPE *) realloc(G_hat, Mesh*Mesh*Mesh*sizeof(G_hat));  

#ifdef SUBSTRACT_SF
  SelfForce_x = (FLOAT_TYPE *) realloc(SelfForce_x, Teilchenzahl*sizeof(FLOAT_TYPE));
  SelfForce_y = (FLOAT_TYPE *) realloc(SelfForce_y, Teilchenzahl*sizeof(FLOAT_TYPE));
 SelfForce_z = (FLOAT_TYPE *) realloc(SelfForce_z, Teilchenzahl*sizeof(FLOAT_TYPE));
#endif

  forward_plan = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)Qmesh, (fftw_complex *)Qmesh, FFTW_FORWARD, FFTW_ESTIMATE);
  backward_plan = fftw_plan_dft_3d(Mesh, Mesh, Mesh, (fftw_complex *)Qmesh, (fftw_complex *)Qmesh, FFTW_BACKWARD, FFTW_ESTIMATE);
  
}

