#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "p3m.h"
#include "interpol.h"

// Methods

#include "p3m-ik-i.h"
#include "p3m-ik.h"
#include "p3m-ad-i.h"
#include "ewald.h"

/* Name der Datei, in der die zu bearbeitende Konfiguration steht: */
// #define EINLESDAT "C100_10" 
 #define EINLESDAT "C4_10" 
/* Name der Datei, in der die exakten Energien und Kraefte stehen: */
// #define EXAKDAT "C100_10.exakt" 
 #define EXAKDAT "C4_10.exakt" 

// #define WRITE_FORCES

// #define FORCE_DEBUG

FLOAT_TYPE *Fx_exa, *Fy_exa, *Fz_exa;
FLOAT_TYPE *Fx, *Fy, *Fz;
FLOAT_TYPE *Fx_R, *Fy_R, *Fz_R;
FLOAT_TYPE *Fx_D, *Fy_D, *Fz_D;
int Teilchenzahl;
int kmax;
FLOAT_TYPE Temp, Bjerrum;
FLOAT_TYPE alpha;
FLOAT_TYPE rcut;
FLOAT_TYPE beta;

void identity(void) {
  return;
}

void (*Influence_function_berechnen)(FLOAT_TYPE);
void (*P3M)(FLOAT_TYPE,int);
void (*Init)(int);


void Realteil(FLOAT_TYPE alpha)
{
  /* Zwei Teilchennummern: */
  int t1,t2;
  /* Minimum-Image-Abstand: */
  FLOAT_TYPE dx,dy,dz,r;
  /* Staerke der elegktrostatischen Kraefte */
  FLOAT_TYPE fak;
  /* Zur Approximation der Fehlerfunktion */
  FLOAT_TYPE t,ar,erfc_teil;

  const FLOAT_TYPE wupi = 1.77245385090551602729816748334;

  /* Zur Approximation der komplementaeren Fehlerfunktion benoetigte
     Konstanten. Die Approximation stammt aus Abramowitz/Stegun:
     Handbook of Mathematical Functions, Dover (9. ed.), Kapitel 7. */
  const FLOAT_TYPE a1 =  0.254829592;
  const FLOAT_TYPE a2 = -0.284496736;
  const FLOAT_TYPE a3 =  1.421413741;
  const FLOAT_TYPE a4 = -1.453152027;
  const FLOAT_TYPE a5 =  1.061405429;
  const FLOAT_TYPE  p =  0.3275911;

#pragma omp paralell for privat(dx,dy,dz, r, ar, erfc_teil, fak)
  for (t1=0; t1<Teilchenzahl-1; t1++)   /* Quick and Dirty N^2 */
    for (t2=t1+1; t2<Teilchenzahl; t2++)
      {
	dx = xS[t1] - xS[t2]; dx -= round(dx*Leni)*Len;
	dy = yS[t1] - yS[t2]; dy -= round(dy*Leni)*Len; 
	dz = zS[t1] - zS[t2]; dz -= round(dz*Leni)*Len;
	r = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
	if (r<=rcut)
	  {
	    ar= alpha*r;
	    //	    t = 1.0 / (1.0 + p*ar);
	    //      erfc_teil = t*(a1+t*(a2+t*(a3+t*(a4+t*a5))));
            erfc_teil = erfc(ar);
	    fak = Q[t1]*Q[t2]*
	      (erfc_teil/r+(2*alpha/wupi)*exp(-ar*ar))/SQR(r);
	    
	    Fx_R[t1] += fak*dx;
	    Fy_R[t1] += fak*dy;
	    Fz_R[t1] += fak*dz;
	    Fx_R[t2] -= fak*dx;
	    Fy_R[t2] -= fak*dy;
	    Fz_R[t2] -= fak*dz;
	  }	 
      }
}


void Dipol(FLOAT_TYPE alpha)
{
  FLOAT_TYPE SummeX,SummeY,SummeZ;
  FLOAT_TYPE VorFak;
  FLOAT_TYPE d;
  int i;

  SummeX=SummeY=SummeZ=0.0;
  VorFak=4.0*PI/(3.0*Len*Len*Len);

  for (i=0; i<Teilchenzahl; ++i)
    if (Q[i]!=0)
      {
	SummeX += Q[i]*xS[i];
	SummeY += Q[i]*yS[i];
	SummeZ += Q[i]*zS[i];
      }

 //   Coulomb-Energie: 

  SummeX *= (Bjerrum*Temp*VorFak);
  SummeY *= (Bjerrum*Temp*VorFak);  
  SummeZ *= (Bjerrum*Temp*VorFak);

 // Kraefte: 
  for (i=0; i<Teilchenzahl; ++i)
    if (Q[i]!=0)
      {
	Fx_D[i] -= Q[i]*SummeX;
	Fy_D[i] -= Q[i]*SummeY;
	Fz_D[i] -= Q[i]*SummeZ;
      }
}

void Elstat_berechnen(FLOAT_TYPE alpha)
{
  /* 
     Zuerst werden die Kraefte und Energien auf Null 
     gesetzt, dann werden die Beitraege zur Elektrostatik aus
     Orts-und Impulsraum sowie die Dipolkorrektur aufaddiert. (P3M)
  */

  int i;
  for (i=0; i<Teilchenzahl; i++) 
    {
      Fx[i]   = Fy[i]   = Fz[i]   = 0.0;
      Fx_R[i] = Fy_R[i] = Fz_R[i] = 0.0;
      Fx_K[i] = Fy_K[i] = Fz_K[i] = 0.0;
      Fx_D[i] = Fy_D[i] = Fz_D[i] = 0.0;
    }
    

  
  E_Coulomb_Dipol = E_Coulomb_Self = E_Coulomb_Impuls_Summe = E_Coulomb_Real_Summe = 0.0;
  
  Realteil(alpha);

  //  Dipol(alpha);
  
  P3M(alpha, Teilchenzahl);
  
  for (i=0; i<Teilchenzahl; i++) 
    {
      Fx[i] += (Fx_R[i]+Fx_K[i]);
      Fy[i] += (Fy_R[i]+Fy_K[i]);
      Fz[i] += (Fz_R[i]+Fz_K[i]);
#ifdef FORCE_DEBUG
      printf("Particle %d Total Force (%g %g %g) [R(%g %g %g) K(%g %g %g)]\n", i, Fx[i], Fy[i], Fz[i], Fx_R[i], Fy_R[i], Fz_R[i], Fx_K[i], Fy_K[i], Fz_K[i]);
#endif
    }
}


void Differenzenoperator_berechnen(void)
{
  /* 
     Die Routine berechnet den fourieretransformierten 
     Differentialoperator auf der Ebene der n, nicht der k,
     d.h. der Faktor  i*2*PI/L fehlt hier!
  */
  
  int    i;
  FLOAT_TYPE dMesh=(FLOAT_TYPE)Mesh;
  FLOAT_TYPE dn;

  Dn = (FLOAT_TYPE *) realloc(Dn, Mesh*sizeof(FLOAT_TYPE));  

  fprintf(stderr,"Fouriertransformierten DIFFERENTIAL-Operator vorberechnen...");
  
  for (i=0; i<Mesh; i++)
    {
      dn    = (FLOAT_TYPE)i; 
      dn   -= round(dn/dMesh)*dMesh;
      Dn[i] = dn;
    }

  Dn[Mesh/2] = 0.0;  

  fprintf(stderr,"\n");
}

void nshift_ausrechnen(void)
{
  /* Verschiebt die Meshpunkte um Mesh/2 */
  
  int    i;
  FLOAT_TYPE dMesh=(FLOAT_TYPE)Mesh;

  fprintf(stderr,"Mesh-Verschiebung vorberechnen...");

  nshift = (FLOAT_TYPE *) realloc(nshift, Mesh*sizeof(FLOAT_TYPE));  

  for (i=0; i<Mesh; i++) nshift[i] = i - round(i/dMesh)*dMesh; 
  
  fprintf(stderr,"\n");
}

void Exakte_Werte_einlesen(char *filename)
{
  /* Liest die exakten Werte fuer Energien und Kraefte der
     einzelnen Teilchen ein. Dateiname muss im #define EXAKTDAT
     angegeben werden. */
  
  FILE *fp;
  int i;
  FLOAT_TYPE E_Coulomb;
  
  fp=fopen(filename, "r");
  
  if((fp == NULL) || feof(fp)) {
    fprintf(stderr, "Could not open '%s' for reading.\n", filename);
  }

  Fx_exa = (FLOAT_TYPE *) realloc(Fx_exa, Teilchenzahl*sizeof(FLOAT_TYPE));
  Fy_exa = (FLOAT_TYPE *) realloc(Fy_exa, Teilchenzahl*sizeof(FLOAT_TYPE));
  Fz_exa = (FLOAT_TYPE *) realloc(Fz_exa, Teilchenzahl*sizeof(FLOAT_TYPE));

  for (i=0; i<Teilchenzahl; ++i)
    fscanf(fp,"%lf\t%lf\t%lf\t%lf\n",
	   &E_Coulomb,
	   &Fx_exa[i],
	   &Fy_exa[i],
	   &Fz_exa[i]); 
  
  fclose(fp);
}

void Daten_einlesen(char *filename)
{
  /* Oeffnet die Datei "PMETest.dat" zum LESEN. Dort muessen
     die systemrelevanten Daten sowie Orte und Ladungen der
     Teilchen stehen. */
  
  FILE *fp;
  int i;
  
  fp=fopen(filename, "r");

  if((fp == NULL) || feof(fp)) {
    fprintf(stderr, "Could not open '%s' for reading.\n", filename);
  }
      

  fscanf(fp,"# Teilchenzahl: %d\n",&Teilchenzahl);
  fscanf(fp,"# Len: %lf\n",&Len);
  fscanf(fp,"# Mesh: %d\n",&Mesh);
  fscanf(fp,"# kmax: %d\n",&kmax);
  fscanf(fp,"# alpha: %lf\n",&alpha);
  fscanf(fp,"# beta: %lf\n",&beta);
  fscanf(fp,"# ip: %d\n",&ip);
  fscanf(fp,"# rcut: %lf\n",&rcut);
  fscanf(fp,"# Temp: %lf\n",&Temp);
  fscanf(fp,"# Bjerrum: %lf\n",&Bjerrum);


  fprintf(stderr,"# Teilchenzahl: %d\n", Teilchenzahl);
  fprintf(stderr,"# Len:          %lf\n",Len);
  fprintf(stderr,"# Mesh:         %d\n", Mesh);
  fprintf(stderr,"# kmax:         %d\n", kmax);
  fprintf(stderr,"# alpha:        %lf\n",alpha);
  fprintf(stderr,"# beta:         %lf\n",beta);
  fprintf(stderr,"# ip:           %d\n", ip);
  fprintf(stderr,"# rcut:         %lf\n",rcut);
  fprintf(stderr,"# Temp:         %lf\n",Temp);
  fprintf(stderr,"# Bjerrum:      %lf\n",Bjerrum);

  xS = (FLOAT_TYPE *) realloc(xS, Teilchenzahl*sizeof(FLOAT_TYPE));
  yS = (FLOAT_TYPE *) realloc(yS, Teilchenzahl*sizeof(FLOAT_TYPE));
  zS = (FLOAT_TYPE *) realloc(zS, Teilchenzahl*sizeof(FLOAT_TYPE));
   Q = (FLOAT_TYPE *) realloc( Q, Teilchenzahl*sizeof(FLOAT_TYPE));

  Fx = (FLOAT_TYPE *) realloc(Fx, Teilchenzahl*sizeof(FLOAT_TYPE));
  Fy = (FLOAT_TYPE *) realloc(Fy, Teilchenzahl*sizeof(FLOAT_TYPE));
  Fz = (FLOAT_TYPE *) realloc(Fz, Teilchenzahl*sizeof(FLOAT_TYPE));

  Fx_R = (FLOAT_TYPE *) realloc(Fx_R, Teilchenzahl*sizeof(FLOAT_TYPE));
  Fy_R = (FLOAT_TYPE *) realloc(Fy_R, Teilchenzahl*sizeof(FLOAT_TYPE));
  Fz_R = (FLOAT_TYPE *) realloc(Fz_R, Teilchenzahl*sizeof(FLOAT_TYPE));

  Fx_K = (FLOAT_TYPE *) realloc(Fx_K, Teilchenzahl*sizeof(FLOAT_TYPE));
  Fy_K = (FLOAT_TYPE *) realloc(Fy_K, Teilchenzahl*sizeof(FLOAT_TYPE));
  Fz_K = (FLOAT_TYPE *) realloc(Fz_K, Teilchenzahl*sizeof(FLOAT_TYPE));

  Fx_D = (FLOAT_TYPE *) realloc(Fx_D, Teilchenzahl*sizeof(FLOAT_TYPE));
  Fy_D = (FLOAT_TYPE *) realloc(Fy_D, Teilchenzahl*sizeof(FLOAT_TYPE));
  Fz_D = (FLOAT_TYPE *) realloc(Fz_D, Teilchenzahl*sizeof(FLOAT_TYPE));

  Leni = 1.0 / Len;
  Q2 = 0.0;
  /* Teilchenkoordinaten und -ladungen: */
  for (i=0; i<Teilchenzahl; i++) {
    fscanf(fp,"%lf\t%lf\t%lf\t%lf\n",&xS[i],&yS[i],&zS[i],&Q[i]);
    Q2 += Q[i]*Q[i];
  }
  fclose(fp);
}

void usage(char *name) {
  fprintf(stderr, "usage: %s <positions> <forces> <alpha_min> <alpha_max> <alpha_step> <method>\n", name);
}

int main(int argc, char **argv)
{
  int    i;
  long   m;
  FLOAT_TYPE d;
  FLOAT_TYPE sx,sy,sz;
  FLOAT_TYPE EC2,DeltaEC2,FC2,DeltaFC2;
  FLOAT_TYPE DeltaE_rel,DeltaE_abs,DeltaF_rel,DeltaF_abs;
  FLOAT_TYPE alphamin,alphamax,alphastep;
  int method;
  char *method_name;  

  FILE* fout;

  if(argc != 7) {
    usage(argv[0]);
    return 128;
  }

  Daten_einlesen(argv[1]);

  nshift_ausrechnen();

  alphamin = atof(argv[3]);
  alphamax = atof(argv[4]);
  alphastep = atof(argv[5]);

  method = atoi(argv[6]);

  Interpolationspolynom_berechnen(ip);    /* Hockney/Eastwood */

  /*
  for (i=0; i<Teilchenzahl; i++)
    {
      xS[i]-=floor(xS[i]*Leni)*Len;
      yS[i]-=floor(yS[i]*Leni)*Len;
      zS[i]-=floor(zS[i]*Leni)*Len;
    }
  */
  Differenzenoperator_berechnen();    /* kontinuierlich, d.h. DIFFERENTIAL! */
  


  #ifdef WRITE_FORCES
  alpha = Ewald_compute_optimal_alpha(rcut, Teilchenzahl);
  printf("Optimal alpha is %lf (error: %e)\n", alpha, Ewald_estimate_error(alpha, rcut, Teilchenzahl));
  Ewald_init(Teilchenzahl);
  Ewald_compute_influence_function(alpha);
  P3M = &Ewald_k_space;
  Elstat_berechnen(alpha);
  fout = fopen(argv[2], "w");
  printf("Writing forces to %s\n", argv[2]);
  for(i=0;i<Teilchenzahl;i++) {
    fprintf(fout, "%d %.22e %.22e %.22e\n", i, Fx[i], Fy[i], Fz[i]);
  }
  fclose(fout);
  #endif

  Exakte_Werte_einlesen(argv[2]);


  switch(method) {
    case 0:
      method_name = "P3M with ik-differentiation, not interlaced";
      P3M = &P3M_ik;
      Influence_function_berechnen = &Influence_function_berechnen_ik;
      Init = &Init_ik;
      break;
    case 1:
      method_name = "P3M with ik-differentiation, interlaced";
      P3M = &P3M_ik_interlaced;
      Influence_function_berechnen = &Influence_function_berechnen_ik_interlaced;
      Init = &Init_interlaced_ik;
      break;
    case 2:
        method_name = "Ewald summation.";
	P3M = &Ewald_k_space;
	Influence_function_berechnen = &Ewald_compute_influence_function;
	Init = &Ewald_init;
        break;
    case 3:
      method_name = "P3M with analytical diff, interlaced.";
      P3M = &P3M_ad_interlaced;
      Influence_function_berechnen = &Influence_function_ad_interlaced;
      Init = &Init_ad_interlaced;
      break;
    default:
      fprintf(stderr, "Method %d not know.\n", method);
      return 127;
  }

  printf("Using %s.\n", method_name);

  //VB
  fout = fopen("out.dat","w");

  Init(Teilchenzahl);

  printf("# alpha\tDeltaF_abs\tDeltaF_rel\n");
  for (alpha=alphamin; alpha<=alphamax; alpha+=alphastep)
    { 
      Influence_function_berechnen(alpha);  /* Hockney/Eastwood */
      
      Elstat_berechnen(alpha); /* Hockney/Eastwood */
      
      /* ACHTUNG: 
	 Der Dipolanteil ist in der Energie NICHT dabei! 
	 Ausserdem wird die Energie zuerst summiert, dann wird 
	 die Differenz zur exakten Gesamtenergie gebildet und dann
	 wird durch die Teilchenzahl dividiert.
      */
      EC2 = FC2 = DeltaFC2 = 0.0;
      for (i=0; i<Teilchenzahl; i++)
	{
#ifdef FORCE_DEBUG
	  printf("%d (%e %e %e), (%e %e %e)\n", i, Fx_exa[i], Fy_exa[i], Fz_exa[i], Fx[i], Fy[i], Fz[i]); 
#endif
	  FC2 += SQR(Fx_exa[i]) + SQR(Fy_exa[i]) + SQR(Fz_exa[i]);
	  DeltaFC2 += SQR(Fx[i]-Fx_exa[i])+SQR(Fy[i]-Fy_exa[i])+SQR(Fz[i]-Fz_exa[i]);
	}
      
      DeltaF_abs = sqrt(DeltaFC2/(FLOAT_TYPE)Teilchenzahl);
      DeltaF_rel = DeltaF_abs/sqrt(FC2);
      
      /* AUSGABE:
	 1. Spalte: alpha
	 2. Spalte: absoluter Fehler in der Kraft
	 3. Spalte: relativer Fehler in der Kraft
      */
      if(method == 2)
        printf("% lf\t% e\t% e\t (est: %e)\n", alpha,DeltaF_abs,DeltaF_rel, Ewald_estimate_error(alpha, rcut, Teilchenzahl));
      else
	printf("% lf\t% e\t% e\n", alpha,DeltaF_abs,DeltaF_rel);
      fflush(stdout);
      fprintf(fout,"% lf\t% e\t% e\n",alpha,DeltaF_abs, DeltaF_rel);

    }
    fclose(fout);

}

