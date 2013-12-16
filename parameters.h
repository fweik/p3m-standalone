/**    Copyright (C) 2011,2012,2013 Florian Weik <fweik@icp.uni-stuttgart.de>

       This program is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

       This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this program.  If not, see <http://www.gnu.org/licenses/>. **/

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
