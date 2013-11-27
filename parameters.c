#include "parameters.h"

void print_parameter(cmd_parameter_t *it) {
  if(it->is_set) 
      switch(it->type) {
        case ARG_TYPE_INT:
	  printf("'%s' = %d, ", it->parameter, *(it->value.i));
	  break;
        case ARG_TYPE_FLOAT:
	  printf("'%s' = %lf, ", it->parameter, FLOAT_CAST *(it->value.f));
	  break;
        case ARG_TYPE_STRING:
	  printf("'%s' = '%s', ", it->parameter, *(it->value.c));
	  break;
        case ARG_TYPE_NONE:
          printf("'%s' is set, ", it->parameter);
      }      
  else
    printf("'%s' = not set, ", it->parameter);
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
      print_parameter(params.required[i]);
    }
  for(i=0;i<params.n_opt;i++)
    if(params.optional[i]->is_set)
      print_parameter(params.optional[i]);

  printf("\n");
}
