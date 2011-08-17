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

// Real space part

#include "realpart.h"

// Dipol correction

#include "dipol.h"

// #define WRITE_FORCES

// #define FORCE_DEBUG
// #define CA_DEBUG

void Elstat_berechnen(system_t *, p3m_parameters_t *, const method_t *);

void Elstat_berechnen(system_t *s, p3m_parameters_t *p, const method_t *m)
{
  /* 
     Zuerst werden die Kraefte und Energien auf Null 
     gesetzt, dann werden die Beitraege zur Elektrostatik aus
     Orts-und Impulsraum sowie die Dipolkorrektur aufaddiert. (P3M)
  */

  int i, j;
  
  for (i=0; i<3; i++) 
    {
      memset(s->f.fields[i], 0, s->nparticles*sizeof(FLOAT_TYPE));
      memset(s->f_k.fields[i], 0, s->nparticles*sizeof(FLOAT_TYPE));
      memset(s->f_r.fields[i], 0, s->nparticles*sizeof(FLOAT_TYPE));
    }
  
  Realteil(s, p);

  //  Dipol(s, p);
  
  m->Kspace_force(s, p);

  for(j=0; j < 3; j++) {
    #pragma omp parallel for
    for (i=0; i<s->nparticles; i++) 	    
      {
	s->f.fields[j][i] += s->f_k.fields[j][i] + s->f_r.fields[j][i];
      }
  }
}

void usage(char *name) {
  fprintf(stderr, "usage: %s <positions> <forces> <alpha_min> <alpha_max> <alpha_step> <method>\n", name);
}

void calc_reference_forces(system_t *s, p3m_parameters_t *p) {
  int i,j;
  p3m_parameters_t op = *p;
  
  op.alpha = Ewald_compute_optimal_alpha(op.rcut, s->nparticles);
  
  fprintf(stderr, "Optimal alpha is %lf (error: %e)\n", op.alpha, Ewald_estimate_error(op.alpha, op.rcut, s->nparticles));
  
  Ewald_init(s->nparticles);
  Ewald_compute_influence_function(op.alpha);

  Elstat_berechnen(s, &op, &method_ewald);

  for(j=0;j<3;j++)
    for(i=0;i<s->nparticles;i++) {
      s->reference.f.fields[j][i] = s->f.fields[j][i];
      s->reference.f_k.fields[j][i] = s->f_k.fields[j][i];
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

  if((method.Init == NULL) || (method.Influence_function == NULL) || (method.Kspace_force == NULL)) {
    fprintf(stderr,"Internal error: Method '%s' (%d) is not correctly defined. Aborting.\n", method.method_name, emthod.method_id);
    exit(-1);
  }

  fprintf(stderr, "Using %s.\n", method.method_name);

  fout = fopen("out.dat","w");

  method.Init(&system, &parameters);

  printf("# %8s\t%8s\t%8s\t%8s\t%8s\n", "alpha", "DeltaF", "Estimate", "R-Error", "K-Error");
  for (parameters.alpha=alphamin; parameters.alpha<=alphamax; parameters.alpha+=alphastep)
    { 
      method.Influence_function(&system, &parameters);  /* Hockney/Eastwood */
  
      Elstat_berechnen(&system, &parameters, &method); /* Hockney/Eastwood */
  
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

