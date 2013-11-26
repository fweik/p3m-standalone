#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "types.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

enum { ARG_TYPE_INT, ARG_TYPE_FLOAT, ARG_TYPE_STRING, ARG_TYPE_NONE };
enum { ARG_OPTIONAL, ARG_REQUIRED };

typedef struct {
  char *parameter;
  int type;
  unsigned int is_set : 1;
  union {
    int *i;
    FLOAT_TYPE *f;
    char **c;
  } value;    
} cmd_parameter_t;

typedef struct {
  cmd_parameter_t **required;
  int n_req;
  cmd_parameter_t **optional; 
  int n_opt;
} cmd_parameters_t;

// functions for command line parameter handling
void add_param( char *name, int type, int reqopt, void *value, cmd_parameters_t *params);
int param_isset(char *s, cmd_parameters_t params);
void parse_parameters( int, char **, cmd_parameters_t);

#endif
