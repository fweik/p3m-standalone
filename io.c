#include <stdio.h>
#include <stdlib.h>

#include "p3m.h"

void Exakte_Werte_einlesen(char *filename, int number_of_particles)
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

  Fx_exa = (FLOAT_TYPE *) realloc(Fx_exa, number_of_particles*sizeof(FLOAT_TYPE));
  Fy_exa = (FLOAT_TYPE *) realloc(Fy_exa, number_of_particles*sizeof(FLOAT_TYPE));
  Fz_exa = (FLOAT_TYPE *) realloc(Fz_exa, number_of_particles*sizeof(FLOAT_TYPE));



  for (i=0; i<number_of_particles; ++i)
    fscanf(fp,"%lf\t%lf\t%lf\t%lf\n",
	   &E_Coulomb,
	   &Fx_exa[i],
	   &Fy_exa[i],
	   &Fz_exa[i]); 
  
  fclose(fp);
}

void init_arrays(int Teilchenzahl) {

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

  if((Mesh & (Mesh - 1)) != 0) {
    fprintf(stderr, "Meshsize must be power of 2!.");
    exit(128);
  }

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

  cao = ip + 1;
  cao3 = cao*cao*cao;

  init_arrays(Teilchenzahl);

  Leni = 1.0 / Len;
  Q2 = 0.0;
  /* Teilchenkoordinaten und -ladungen: */
  for (i=0; i<Teilchenzahl; i++) {
    fscanf(fp,"%lf\t%lf\t%lf\t%lf\n",&xS[i],&yS[i],&zS[i],&Q[i]);
    Q2 += Q[i]*Q[i];
  }
  fclose(fp);
}

