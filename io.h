#ifndef IO_H
#define IO_H

#include "types.h"

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

void Read_exact_forces(system_t *s, char *);
system_t *Read_system(parameters_t *, char *);
void Write_exact_forces(system_t *, char *);
void Write_system(system_t *, char *);
void Write_system_cuda( system_t *s, parameters_t *p, char *filename);

void write_vtf(char *filename, system_t *s);

#endif
