#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "p3m.h"
#include "interpol.h"

// Methods

#include "p3m-ik-i.h"
#include "p3m-ik.h"
#include "p3m-ad.h"
#include "p3m-ad-i.h"
#include "ewald.h"

// Error estimates for p3m

#include "p3m-error.h"

// Utils and IO

#include "io.h"

/* Name der Datei, in der die zu bearbeitende Konfiguration steht: */
// #define EINLESDAT "C100_10" 
 #define EINLESDAT "C4_10" 
/* Name der Datei, in der die exakten Energien und Kraefte stehen: */
// #define EXAKDAT "C100_10.exakt" 
 #define EXAKDAT "C4_10.exakt" 

// #define WRITE_FORCES

 #define FORCE_DEBUG
// #define CA_DEBUG

#define r_ind(A,B,C) ((A)*Mesh*Mesh + (B)*Mesh + (C))
#define c_ind(A,B,C) (2*Mesh*Mesh*(A)+2*Mesh*(B)+2*(C))

void print_ghat() {
  int i;
  if(G_hat == NULL)
    return;
  for(i=0;i<(Mesh*Mesh*Mesh);i++)
    printf("G_hat %e\n", G_hat[i]);
}

void (*Influence_function_berechnen)(FLOAT_TYPE);
void (*kspace_force)(FLOAT_TYPE,int);
void (*Init)(int);
double (*error)(double, int*, int, int, double, double, double, double *) = NULL;


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

#pragma omp paralell for privat(dx,dy,dz, r, ar, erfc_teil, fak, t2)
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
	    //t = 1.0 / (1.0 + p*ar);
	    //erfc_teil = t*(a1+t*(a2+t*(a3+t*(a4+t*a5))));
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
  
  kspace_force(alpha, Teilchenzahl);
  
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

void usage(char *name) {
  fprintf(stderr, "usage: %s <positions> <forces> <alpha_min> <alpha_max> <alpha_step> <method>\n", name);
}

void check_g(void) {
  int i,j,k;
  FILE *f = fopen("p3m-wall-gforce.dat", "r");
  double g_now = 0.0, rms_g = 0.0;

  printf("g handle %p\n", f);

  for(i=0;i<64;i++) 
    for(j=0;j<64;j++)
      for(k=0;k<64;k++) {
	fscanf( f, "%lf", &g_now);
	//	printf("g_diff %e %e %e\n", G_hat[r_ind(k,i,j)], g_now, G_hat[r_ind(k,i,j)] - g_now);
	rms_g += SQR(G_hat[r_ind(k,i,j)] - g_now);
  }
  printf("g_hat rms %e\n", sqrt(rms_g)/(64*64*64));

}

void check_k(void) {
  FILE *f = fopen("p3m-wall-k.dat", "r");

  double rms_k = 0.0, rms_r =0.0;
    int i;
  FLOAT_TYPE E_Coulomb;
  
  Fxk_exa = (FLOAT_TYPE *) realloc(Fxk_exa, Teilchenzahl*sizeof(FLOAT_TYPE));
  Fyk_exa = (FLOAT_TYPE *) realloc(Fyk_exa, Teilchenzahl*sizeof(FLOAT_TYPE));
  Fzk_exa = (FLOAT_TYPE *) realloc(Fzk_exa, Teilchenzahl*sizeof(FLOAT_TYPE));



  for (i=0; i<Teilchenzahl; ++i) {
    fscanf(f,"%lf\t%lf\t%lf\t%lf\n",
	   &E_Coulomb,
	   &Fxk_exa[i],
	   &Fyk_exa[i],
	   &Fzk_exa[i]); 

    rms_k += SQR(Fxk_exa[i] - Fx_K[i]);
    rms_k += SQR(Fyk_exa[i] - Fy_K[i]);
    rms_k += SQR(Fzk_exa[i] - Fz_K[i]);

    rms_r += SQR(Fx_exa[i] - Fxk_exa[i] - Fx_R[i]);
    rms_r += SQR(Fy_exa[i] - Fyk_exa[i] - Fy_R[i]);
    rms_r += SQR(Fz_exa[i] - Fzk_exa[i] - Fz_R[i]);
  }
  fclose(f);

  printf("kspace deviation %e\n", sqrt(rms_k)/Teilchenzahl);
  printf("rspace deviation %e\n", sqrt(rms_r)/Teilchenzahl);
}

int main(int argc, char **argv)
{
  int    i;
  long   m;
  FLOAT_TYPE d;
  FLOAT_TYPE sx,sy,sz;
  FLOAT_TYPE EC2,DeltaEC2,FC2,DeltaFC2;
  FLOAT_TYPE rms_x=0.0, rms_y=0.0, rms_z=0.0;
  FLOAT_TYPE DeltaE_rel,DeltaE_abs,DeltaF_rel,DeltaF_abs;
  FLOAT_TYPE alphamin,alphamax,alphastep;
  FLOAT_TYPE rms_k = 0.0, rms_r = 0.0;
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
  kspace_force = &Ewald_k_space;
  Elstat_berechnen(alpha);
  fout = fopen(argv[2], "w");
  printf("Writing forces to %s\n", argv[2]);
  for(i=0;i<Teilchenzahl;i++) {
    fprintf(fout, "%d %.22e %.22e %.22e %.22e %.22e %.22e\n", 
	    i, Fx[i], Fy[i], Fz[i], Fx_K[i], Fy_K[i], Fz_K[i]);
  }
  fclose(fout);
  #endif

  Exakte_Werte_einlesen(argv[2], Teilchenzahl, &Fx_exa, &Fy_exa, &Fz_exa);

  //  printf("Teilchenzahl %d, F_exa ( %p %p %p )\n", Teilchenzahl, Fx_exa, Fy_exa, Fz_exa);

  switch(method) {
    case 0:
      method_name = "P3M with ik-differentiation, not interlaced";
      kspace_force = &P3M_ik;
      Influence_function_berechnen = &Influence_function_berechnen_ik;
      Init = &Init_ik;
      error = &p3m_error_ik;
      break;
    case 1:
      method_name = "P3M with ik-differentiation, interlaced";
      kspace_force = &P3M_ik_interlaced;
      Influence_function_berechnen = &Influence_function_berechnen_ik_interlaced;
      Init = &Init_interlaced_ik;
      break;
    case 2:
      method_name = "Ewald summation.";
      kspace_force = &Ewald_k_space;
      Influence_function_berechnen = &Ewald_compute_influence_function;
      Init = &Ewald_init;
      error = &Ewald_error_wrapper;
      break;
    case 3:
      method_name = "P3M with analytical diff, noninterlaced.";
      kspace_force = &P3M_ad;
      Influence_function_berechnen = &Influence_function_berechnen_ad;
      Init = &Init_ad;
      error = &p3m_error_ad;
      break;
    case 4:
      method_name = "P3M with analytical diff, interlaced.";
      kspace_force = &P3M_ad_interlaced;
      Influence_function_berechnen = &Influence_function_berechnen_ad_interlaced;
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

  printf("# %8s\t%8s\t%8s\t%8s\t%8s\n", "alpha", "DeltaF", "Estimate", "R-Error", "K-Error");
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
      rms_k = rms_r = 0.0;
      for (i=0; i<Teilchenzahl; i++)
	{
	  #ifdef FORCE_DEBUG
          printf("Particle %d Total Force (%g %g %g) [R(%g %g %g) K(%g %g %g)] reference (%g %g %g)\n", i, Fx[i], Fy[i], Fz[i], Fx_R[i], Fy_R[i], Fz_R[i], Fx_K[i], Fy_K[i], Fz_K[i], Fx_exa[i], Fy_exa[i], Fz_exa[i]);
          rms_x += SQR(Fx_exa[i]);
          rms_y += SQR(Fy_exa[i]);
          rms_z += SQR(Fz_exa[i]);
#endif // FORCE_DEBUG
	  rms_k += SQR(Fxk_exa[i] - Fx_K[i]);
	  rms_k += SQR(Fyk_exa[i] - Fy_K[i]);
	  rms_k += SQR(Fzk_exa[i] - Fz_K[i]);
	  
	  rms_r += SQR(Fx_exa[i] - Fxk_exa[i] - Fx_R[i]);
	  rms_r += SQR(Fy_exa[i] - Fyk_exa[i] - Fy_R[i]);
	  rms_r += SQR(Fz_exa[i] - Fzk_exa[i] - Fz_R[i]);

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
      if(error != NULL) {
	int mesh[3] = {Mesh, Mesh, Mesh};
        double box[3] = {Len, Len, Len};
        double estimate =  error(1.0, mesh, cao, Teilchenzahl, Q2, alpha*Len, rcut*Leni, box);
        printf("%8lf\t%8e\t%8e\t %8e %8e\n", alpha,DeltaF_abs, estimate, 
	       sqrt(rms_r)/Teilchenzahl, sqrt(rms_k)/Teilchenzahl);
	fprintf(fout,"% lf\t% e\t% e\t% e\n",alpha,DeltaF_abs, DeltaF_rel, estimate);
      }
      else {
	printf("%8lf\t%8e\t na\t%8e\t%8e\n", alpha,DeltaF_abs, sqrt(rms_r)/Teilchenzahl, sqrt(rms_k)/Teilchenzahl);
	fprintf(fout,"% lf\t% e\t% e\t na\n",alpha,DeltaF_abs, DeltaF_rel);
      }
#ifdef FORCE_DEBUG
      fprintf(stderr, "%lf rms %e %e %e\n", alpha, sqrt(rms_x), sqrt(rms_y), sqrt(rms_z));
#endif
      fflush(stdout);
      fflush(fout);
    }
    fclose(fout);

}

