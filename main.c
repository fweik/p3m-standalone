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

// #define WRITE_FORCES

// #define FORCE_DEBUG
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

/*
void (*Influence_function_berechnen)(FLOAT_TYPE);
void (*kspace_force)(FLOAT_TYPE,int);
void (*Init)(int);
double (*error)(double, int*, int, int, double, double, double, double *) = NULL;
*/

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

#pragma omp parallel for private(dx,dy,dz, r, ar, erfc_teil, fak, t2)
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

void usage(char *name) {
  fprintf(stderr, "usage: %s <positions> <forces> <alpha_min> <alpha_max> <alpha_step> <method>\n", name);
}

/*
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
*/
void calc_reference_forces(char *forces_file) {
  int i;

  alpha = Ewald_compute_optimal_alpha(rcut, Teilchenzahl);
  fprintf(stderr, "Optimal alpha is %lf (error: %e)\n", alpha, Ewald_estimate_error(alpha, rcut, Teilchenzahl));
  Ewald_init(Teilchenzahl);
  Ewald_compute_influence_function(alpha);
  kspace_force = &Ewald_k_space;
  Elstat_berechnen(alpha);

  for(i=0;i<Teilchenzahl;i++) {
    Fx_exa[i] = Fx[i];
    Fy_exa[i] = Fy[i];
    Fz_exa[i] = Fz[i];
    Fxk_exa[i] = Fx_K[i];
    Fyk_exa[i] = Fy_K[i];
    Fzk_exa[i] = Fz_K[i];
  }
}

int main(int argc, char **argv)
{
  int    i,j;
  int methodnr;
  FLOAT_TYPE EC2,FC2,DeltaFC2;

  FLOAT_TYPE rms_v[3];

  FLOAT_TYPE DeltaF_rel,DeltaF_abs;
  FLOAT_TYPE alphamin,alphamax,alphastep;
  FLOAT_TYPE rms_k = 0.0, rms_r = 0.0;

  FILE* fout;

  system_t system;
  method_t method;
  p3m_parameters_t parameters;

  if(argc != 7) {
    usage(argv[0]);
    return 128;
  }

  // Inits the system and reads particle data and parameters from file.
  Daten_einlesen(&system, &parameters, argv[1]);

  nshift_ausrechnen();

  alphamin = atof(argv[3]);
  alphamax = atof(argv[4]);
  alphastep = atof(argv[5]);

  methodnr = atoi(argv[6]);

  Interpolationspolynom_berechnen(parameters.ip);    /* Hockney/Eastwood */

  #ifdef WRITE_FORCES
  calc_reference_forces(argv[2]);
  #endif
  #ifndef WRITE_FORCES
  Exakte_Werte_einlesen(&system, argv[2]);
  #endif

  if(methodnr == method_p3m_ik.method_id) 
    method = method_p3m_ik;
  else if(methodnr == method_p3m_ik_i.method_id) 
    method = method_p3m_ik_i;
  else if(methodnr == method_p3m_ad.method_id) 
    method = method_p3m_ad;
  else if(methodnr == method_p3m_ad_i.method_id) 
    method = method_p3m_ad_i;
  else {
    fprintf(stderr, "Method %d not know.", methodnr);
    exit(126);
  }

  fprintf(stderr, "Using %s.\n", method.method_name);

  fout = fopen("out.dat","w");

  method.Init(&system, &parameters);

  printf("# %8s\t%8s\t%8s\t%8s\t%8s\n", "alpha", "DeltaF", "Estimate", "R-Error", "K-Error");
  for (parameters.alpha=alphamin; parameters.alpha<=alphamax; parameters.alpha+=alphastep)
    { 
      method.Influence_function(&system, &parameters);  /* Hockney/Eastwood */
  
      Elstat_berechnen(parameters.alpha); /* Hockney/Eastwood */
  
      EC2 = FC2 = DeltaFC2 = 0.0;
      rms_k = rms_r = 0.0;
      rms_v[0] = rms_v[1] = rms_v[2] = 0.0;
      for (i=0; i<system.nparticles; i++)
	{
	  for(j=0;j<3;j++) {
	    rms_v[j] += SQR(system.reference.f.fields[j][i] - system.f.fields[j][i]);

	    rms_k += SQR(system.reference.f_k.fields[j][i] 
			 - system.f_k.fields[j][i]);
	    rms_r += SQR(system.reference.f.fields[j][i] 
			 - system.reference.f_k.fields[j][i] 
			 - system.f.fields[j][i]);

	    DeltaFC2 += SQR(system.f.fields[j][i] 
			    - system.reference.f.fields[j][i]);
	    
	    FC2 += SQR(system.reference.f.fields[j][i]);
	  }
#ifdef FORCE_DEBUG
          printf("Particle %d Total Force (%g %g %g) [R(%g %g %g) K(%g %g %g)] reference (%g %g %g)\n", i, 
		 Fx[i], Fy[i], Fz[i], Fx_R[i], Fy_R[i], Fz_R[i], Fx_K[i], Fy_K[i], Fz_K[i], Fx_exa[i], Fy_exa[i], Fz_exa[i]);
#endif // FORCE_DEBUG
	}
      
      DeltaF_abs = sqrt(DeltaFC2/(FLOAT_TYPE)system.nparticles);
      DeltaF_rel = DeltaF_abs/sqrt(FC2);
      
      if(method.Error != NULL) {
        double estimate =  method.Error(&system, &parameters);
        printf("%8lf\t%8e\t%8e\t %8e %8e\n", parameters.alpha,DeltaF_abs, estimate, 
	       sqrt(rms_r)/system.nparticles, sqrt(rms_k)/system.nparticles);
	fprintf(fout,"% lf\t% e\t% e\t% e\n",parameters.alpha,DeltaF_abs, DeltaF_rel, estimate);
      }
      else {
	printf("%8lf\t%8e\t na\t%8e\t%8e\n", parameters.alpha,DeltaF_abs, sqrt(rms_r)/system.nparticles, sqrt(rms_k)/system.nparticles);
	fprintf(fout,"% lf\t% e\t% e\t na\n",parameters.alpha,DeltaF_abs, DeltaF_rel);
      }
#ifdef FORCE_DEBUG
      fprintf(stderr, "%lf rms %e %e %e\n", alpha, sqrt(rms_x), sqrt(rms_y), sqrt(rms_z));
#endif
      fflush(stdout);
      fflush(fout);
    }
    fclose(fout);
    
    return 0;
}

