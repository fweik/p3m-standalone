#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "p3m.h"
#include "io.h"

void Exakte_Werte_einlesen(system_t *s, char *filename)
{
    FILE *fp;
    int i;
    FLOAT_TYPE E_Coulomb;

    fp=fopen(filename, "r");

    if ((fp == NULL) || feof(fp)) {
        fprintf(stderr, "Could not open '%s' for reading.\n", filename);
        exit(127);
    }

    for (i=0; i<s->nparticles; i++)
        fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
               &E_Coulomb,
               &s->reference->f->x[i], &s->reference->f->y[i], &s->reference->f->z[i],
               &s->reference->f_k->x[i], &s->reference->f_k->y[i], &s->reference->f_k->z[i]);

    fclose(fp);
}

system_t *Daten_einlesen(parameters_t *p, char *filename)
{
    /* Opens file 'filename' for reanding and reads system parameters,
     particle positions and charges. */

    FILE *fp;
    int i;

    FLOAT_TYPE Temp, Bjerrum, Length;
    int n;

    system_t *s;

    assert(p != NULL);
    assert(filename != NULL);

    fp=fopen(filename, "r");

    if ((fp == NULL) || feof(fp)) {
        fprintf(stderr, "Could not open '%s' for reading.\n", filename);
        exit(127);
    }

    //read system parameters
    fscanf(fp,"# Teilchenzahl: %d\n",&n);
    fscanf(fp,"# Len: %lf\n",&Length);

    //read p3m/ewal parameters
    fscanf(fp,"# Mesh: %d\n",&(p->mesh));
    fscanf(fp,"# alpha: %lf\n",&(p->alpha));
    fscanf(fp,"# ip: %d\n",&(p->ip));
    fscanf(fp,"# rcut: %lf\n",&(p->rcut));
    fscanf(fp,"# Temp: %lf\n",&Temp);
    fscanf(fp,"# Bjerrum: %lf\n",&Bjerrum);

    if ((p->mesh & (p->mesh - 1)) != 0) {
        fprintf(stderr, "Meshsize must be power of 2!.");
        exit(128);
    }

    p->prefactor = Temp*Bjerrum;

    p->cao = p->ip + 1;
    p->cao3 = p->cao*p->cao*p->cao;

    s = Init_system(n);
    s->length = Length;

    s->q2 = 0.0;
    /* Teilchenkoordinaten und -ladungen: */
    for (i=0; i<s->nparticles; i++) {
        fscanf(fp,"%lf\t%lf\t%lf\t%lf\n",&(s->p->x[i]), &(s->p->y[i]), &(s->p->z[i]), &(s->q[i]));
        s->q2 += SQR(s->q[i]);
    }
    fclose(fp);

    return s;
}

void Write_exact_forces(system_t *s, char *forces_file) {
    FILE *fin;
    int i;

    fin = fopen(forces_file, "w");

    if (fin == NULL) {
        fprintf(stderr, "Could not open '%s' for writing!.\n", forces_file);
        exit(127);
    }

    for (i=0;i<s->nparticles;i++) {
        fprintf(fin, "%d %.22e %.22e %.22e %.22e %.22e %.22e\n",
                i, s->reference->f->x[i], s->reference->f->y[i], s->reference->f->z[i],
                s->reference->f_k->x[i], s->reference->f_k->y[i], s->reference->f_k->z[i]);
    }

    fclose(fin);
}
