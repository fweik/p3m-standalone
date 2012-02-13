#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "types.h"
#include "io.h"
#include "common.h"

static void usage(cmd_parameter_t **required, cmd_parameter_t **optinal) {
  
}

int cmd_parameter_t_cmp ( cmd_parameter_t **a, cmd_parameter_t **b ) {
  return strcmp((*a)->parameter, (*b)->parameter);
}

void add_param( char *name, int type, int reqopt, void *value, cmd_parameters_t *params) {
  cmd_parameter_t **list;
  int n;
  cmd_parameter_t *new_param;
  
  if(reqopt == ARG_OPTIONAL) {
    n = params->n_opt;
    params->optional = realloc(params->optional, ++n * sizeof(cmd_parameter_t *));
    list = params->optional;
    params->n_opt = n;
  }
  else if(reqopt == ARG_REQUIRED) {
    n = params->n_req;
    params->required = realloc(params->required, ++n * sizeof(cmd_parameter_t *));
    list = params->required;
    params->n_req = n;
  }

  new_param = malloc(sizeof(cmd_parameter_t));
  list[n-1] = new_param;
  new_param->parameter = name;
  new_param->is_set = 0;
  new_param->type = type;

  switch(type) {
    case ARG_TYPE_INT:
      new_param->value.i = (int *)value;
      break;
    case ARG_TYPE_FLOAT:
      new_param->value.f = (FLOAT_TYPE *)value;
      break;
    case ARG_TYPE_STRING:
      new_param->value.c = (char *)value;
      break;
  }

}

void parse_parameters( int argc, char **argv, cmd_parameters_t params) {
  int i;
  cmd_parameter_t **it, *s;
  cmd_parameter_t search_term;
  s = &search_term;

  printf("Parsing for %d required and %d optional parameters.\n", params.n_req, params.n_opt);
  if((params.n_req == 0) && (params.n_opt == 0)) //nothing to parse
    return;

  //  if(params.n_req > argc) {
  // usage(params.required, params.optional);
  //  exit(23);
  // }

  if(params.n_req > 0)
    qsort(params.required, params.n_req, sizeof(cmd_parameter_t *), (__compar_fn_t)cmd_parameter_t_cmp);
  if(params.n_opt > 0)
    qsort(params.optional, params.n_opt, sizeof(cmd_parameter_t *), (__compar_fn_t)cmd_parameter_t_cmp);
  
  for(i=0;i<params.n_req;i++)
    printf("%s\n", params.required[i]->parameter);

  while(argc > 0) {
    printf("Parsing '%s'.\n", *argv);
    search_term.parameter = *argv;
    it = bsearch(&s, params.required, params.n_req, sizeof(cmd_parameter_t *), (__compar_fn_t)cmd_parameter_t_cmp);

    if(it != NULL) {
      printf("Found '%s', value '%s'\n", (*it)->parameter, *(argv + 1));
      switch((*it)->type) {
        case ARG_TYPE_INT:
	  *((*it)->value.i) = atoi(*(++argv));
	  (*it)->is_set = 1;
	  argc--;
	  break;
        case ARG_TYPE_FLOAT:
	  *((*it)->value.f) = atof(*(++argv));
	  (*it)->is_set = 1;	       
	  printf("Value %lf\n", *((*it)->value.f));
	  break;
        case ARG_TYPE_STRING:
	  (*it)->value.c = *(++argv);
	  (*it)->is_set = 1;
	  argc--;
	  break;
      }      
    }

    argv++;
    argc--;
  }
}

void Read_exact_forces(system_t *s, char *filename)
{
    FILE *fp;
    int i, ret_val;
    FLOAT_TYPE E_Coulomb;

    fp=fopen(filename, "r");

    if ((fp == NULL) || feof(fp)) {
        fprintf(stderr, "Could not open '%s' for reading.\n", filename);
        exit(127);
    }

    for (i=0; i<s->nparticles; i++) {
        ret_val = fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
               &E_Coulomb,
               &s->reference->f->x[i], &s->reference->f->y[i], &s->reference->f->z[i],
               &s->reference->f_k->x[i], &s->reference->f_k->y[i], &s->reference->f_k->z[i]);
	if((ret_val != 7) && (ret_val != 4)) {
          fprintf(stderr, "Error while reading file '%s' (%d)\n", filename, ret_val);
          exit(-1);
        }
    }
    fclose(fp);
}

system_t *Read_system(parameters_t *p, char *filename)
{
    /* Opens file 'filename' for reanding and reads system parameters,
     particle positions and charges. */

    FILE *fp;
    int i;

    FLOAT_TYPE Temp, Bjerrum, Length;
    int n;

    int ret_val = 0;

    system_t *s;

    assert(p != NULL);
    assert(filename != NULL);

    fp=fopen(filename, "r");

    if ((fp == NULL) || feof(fp)) {
        fprintf(stderr, "Could not open '%s' for reading.\n", filename);
        exit(127);
    }

    //read system parameters
    ret_val += fscanf(fp,"# Teilchenzahl: %d\n",&n);
    ret_val += fscanf(fp,"# Len: %lf\n",&Length);

    //read p3m/ewal parameters
    ret_val += fscanf(fp,"# Mesh: %d\n",&(p->mesh));
    ret_val += fscanf(fp,"# alpha: %lf\n",&(p->alpha));
    ret_val += fscanf(fp,"# ip: %d\n",&(p->ip));
    ret_val += fscanf(fp,"# rcut: %lf\n",&(p->rcut));
    ret_val += fscanf(fp,"# Temp: %lf\n",&Temp);
    ret_val += fscanf(fp,"# Bjerrum: %lf\n",&Bjerrum);

    if(ret_val != 8) {
      fprintf(stderr, "Error while reading file '%s'\n", filename);
      exit(-1);
    }
    
    p->prefactor = Temp*Bjerrum;

    p->cao = p->ip + 1;
    p->cao3 = p->cao*p->cao*p->cao;

    s = Init_system(n);
    s->length = Length;

    s->q2 = 0.0;
    /* Teilchenkoordinaten und -ladungen: */
    for (i=0; i<s->nparticles; i++) {
        ret_val = fscanf(fp,"%lf\t%lf\t%lf\t%lf\n",&(s->p->x[i]), &(s->p->y[i]), &(s->p->z[i]), &(s->q[i]));
        if(ret_val != 4) {
          fprintf(stderr, "Error while reading file '%s'\n", filename);
          exit(-1);
        }
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

void Write_system(system_t *s, char *filename)
{
    /* Opens file 'filename' for reanding and reads system parameters,
     particle positions and charges. */

    FILE *fp;
    int i;

    assert(filename != NULL);

    fp=fopen(filename, "w");

    if ((fp == NULL) || feof(fp)) {
        fprintf(stderr, "Could not open '%s' for writing.\n", filename);
        exit(127);
    }

    //read system parameters
    fprintf(fp,"# Teilchenzahl: %d\n", s->nparticles);
    fprintf(fp,"# Len: %lf\n", s->length);

    //read p3m/ewal parameters
    fprintf(fp,"# Mesh: %d\n", 1);
    fprintf(fp,"# alpha: %lf\n", 1.0);
    fprintf(fp,"# ip: %d\n", 5);
    fprintf(fp,"# rcut: %lf\n", 1.0);
    fprintf(fp,"# Temp: %lf\n", 1.0);
    fprintf(fp,"# Bjerrum: %lf\n", 1.0);

    /* Teilchenkoordinaten und -ladungen: */
    for (i=0; i<s->nparticles; i++) {
        fprintf(fp,"%lf\t%lf\t%lf\t%lf\n",(s->p->x[i]), (s->p->y[i]), (s->p->z[i]), (s->q[i]));
    }
    fclose(fp);
}
