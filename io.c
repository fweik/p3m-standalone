#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "types.h"
#include "io.h"
#include "common.h"

static void usage(cmd_parameter_t **required, cmd_parameter_t **optinal) {
  
}

void print_parameter(cmd_parameter_t *it) {
  if(it->is_set) 
      switch(it->type) {
        case ARG_TYPE_INT:
	  printf("'%s' = %d", it->parameter, *(it->value.i));
	  break;
        case ARG_TYPE_FLOAT:
	  printf("'%s' = %lf", it->parameter, FLOAT_CAST *(it->value.f));
	  break;
        case ARG_TYPE_STRING:
	  printf("'%s' = '%s'", it->parameter, *(it->value.c));
	  break;
        case ARG_TYPE_NONE:
          printf("'%s' is set", it->parameter);
      }      
  else
    printf("'%s' = not set", it->parameter);
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
      new_param->value.c = (char **)value;
      break;
  }

}

int param_isset(char *c, cmd_parameters_t params) {
  cmd_parameter_t **it, *s;
  cmd_parameter_t search_term;

  search_term.parameter = c;
  s = &search_term;
  it = bsearch( &s, params.optional, params.n_opt, sizeof(cmd_parameter_t *), (__compar_fn_t)cmd_parameter_t_cmp);

  if(it == NULL)
    return 0;

  return (*it)->is_set;
}

void parse_parameters( int argc, char **argv, cmd_parameters_t params) {
  int i;
  cmd_parameter_t **it, *s;
  cmd_parameter_t search_term;
  s = &search_term;

  printf("Parsing for %d required and %d optional parameters.\n", params.n_req, params.n_opt);
  if((params.n_req == 0) && (params.n_opt == 0)) //nothing to parse
    return;

  if(params.n_req > 0)
    qsort(params.required, params.n_req, sizeof(cmd_parameter_t *), (__compar_fn_t)cmd_parameter_t_cmp);
  if(params.n_opt > 0)
    qsort(params.optional, params.n_opt, sizeof(cmd_parameter_t *), (__compar_fn_t)cmd_parameter_t_cmp);
  
  while(argc > 0) {
    search_term.parameter = *argv;

    it = bsearch(&s, params.required, params.n_req, sizeof(cmd_parameter_t *), (__compar_fn_t)cmd_parameter_t_cmp);
    if(it == NULL)
          it = bsearch(&s, params.optional, params.n_opt, sizeof(cmd_parameter_t *), (__compar_fn_t)cmd_parameter_t_cmp);
    if(it != NULL) {
      switch((*it)->type) {
        case ARG_TYPE_INT:
	  *((*it)->value.i) = atoi(*(++argv));
	  (*it)->is_set = 1;
	  argc--;
	  break;
        case ARG_TYPE_FLOAT:
	  *((*it)->value.f) = atof(*(++argv));
	  (*it)->is_set = 1;	       
	  argc--;
	  break;
        case ARG_TYPE_STRING:
	  *((*it)->value.c) = *(++argv);
	  (*it)->is_set = 1;
	  argc--;
	  break;
        case ARG_TYPE_NONE:
	  (*it)->is_set = 1;
	  break;
      }      
    } else {
      printf("Unknown parameter '%s'\n", *argv);
      exit(23);
    }

    argv++;
    argc--;
  }

  for(i=0;i<params.n_req;i++)
    if(!params.required[i]->is_set) {
      printf("Required parameter '%s' missing.\n", params.required[i]->parameter);
      exit(-1);
    } else {
      print_parameter(params.required[i]); printf(", ");
    }
  for(i=0;i<params.n_opt;i++)
    if(!params.optional[i]->is_set) {
      print_parameter(params.optional[i]); printf(", ");
    }
  printf("\n");
}


void Read_exact_forces(system_t *s, char *filename)
{
    FILE *fp;
    int i, ret_val;
    FLOAT_TYPE E_Coulomb;
    double buf[6];

    fp=fopen(filename, "r");

    if ((fp == NULL) || feof(fp)) {
        fprintf(stderr, "Could not open '%s' for reading.\n", filename);
        exit(127);
    }

    for (i=0; i<s->nparticles; i++) {
        ret_val = fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
			 &E_Coulomb, buf, buf + 1, buf + 2, buf + 3, buf + 4, buf + 5);
	s->reference->f->x[i] = buf[0];
	s->reference->f->y[i] = buf[1];
	s->reference->f->z[i] = buf[2];
	s->reference->f_k->x[i] = buf[3];
	s->reference->f_k->y[i] = buf[4];
	s->reference->f_k->z[i] = buf[5];

			 
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

    double buf[4];
    FLOAT_TYPE Length;
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
    ret_val += fscanf(fp,"# Len: %lf\n",buf);
    Length = buf[0];

    if(ret_val != 2) {
      fprintf(stderr, "Error while reading file '%s'\n", filename);
      exit(-1);
    }
    
    p->cao = p->ip + 1;
    p->cao3 = p->cao*p->cao*p->cao;

    s = Init_system(n);
    s->length = Length;

    s->q2 = 0.0;
    /* Teilchenkoordinaten und -ladungen: */
    for (i=0; i<s->nparticles; i++) {
      ret_val = fscanf(fp,"%lf\t%lf\t%lf\t%lf\n", buf, buf + 1, buf + 2, buf + 3 );
      s->p->x[i] = buf[0];
      s->p->y[i] = buf[1];
      s->p->z[i] = buf[2];
      s->q[i] = buf[3];
	
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
                i, FLOAT_CAST s->reference->f->x[i], FLOAT_CAST s->reference->f->y[i], FLOAT_CAST s->reference->f->z[i],
                FLOAT_CAST s->reference->f_k->x[i], FLOAT_CAST s->reference->f_k->y[i], FLOAT_CAST s->reference->f_k->z[i]);
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
    fprintf(fp,"# Len: %lf\n", FLOAT_CAST  s->length);

    //read p3m/ewal parameters
    fprintf(fp,"# Mesh: %d\n", 1);
    fprintf(fp,"# alpha: %lf\n", 1.0);
    fprintf(fp,"# ip: %d\n", 5);
    fprintf(fp,"# rcut: %lf\n", 1.0);
    fprintf(fp,"# Temp: %lf\n", 1.0);
    fprintf(fp,"# Bjerrum: %lf\n", 1.0);

    /* Teilchenkoordinaten und -ladungen: */
    for (i=0; i<s->nparticles; i++) {
        fprintf(fp,"%lf\t%lf\t%lf\t%lf\n", FLOAT_CAST (s->p->x[i]), FLOAT_CAST (s->p->y[i]), FLOAT_CAST (s->p->z[i]), FLOAT_CAST (s->q[i]));
    }
    fclose(fp);
}
