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

#include "types.h"
#include "io.h"
#include "interpol.h"
#include "charge-assign.h"
#include "p3m-common.h"
#include "common.h"
#include "p3m-ik.c"

int main(int argc, char **argv) {

  cmd_parameters_t params = { NULL, 0, NULL, 0 };
  char *out_file, *pos_file;

  parameters_t parameters;
  system_t *system;
  data_t *data;
  method_t m = { 255, "assign only", METHOD_FLAG_Qmesh | METHOD_FLAG_ca, NULL, NULL, NULL, NULL };

  double *real_data, *F_mesh_real;
  int nm, i;
  int dims[3]; 
  double spacing[3] = { 1.0, 1.0, 1.0 };

  add_param( "mesh", ARG_TYPE_INT, ARG_REQUIRED, &(parameters.mesh), &params );
  add_param( "cao", ARG_TYPE_INT, ARG_REQUIRED, &(parameters.cao), &params );
  add_param( "outfile", ARG_TYPE_STRING, ARG_REQUIRED, &out_file, &params );

  add_param( "positions", ARG_TYPE_STRING, ARG_REQUIRED, &pos_file, &params );

  parse_parameters( argc - 1, argv + 1, params );

  parameters.ip = parameters.cao - 1;
  parameters.cao3 = parameters.cao*parameters.cao*parameters.cao ;

  system = Read_system(&parameters, pos_file);

  data = Init_data(&m, system, &parameters);

  nm = parameters.mesh*parameters.mesh*parameters.mesh;
  real_data = Init_array( nm , sizeof(double));
  F_mesh_real = Init_array( nm, 2*3*sizeof(double));

  dims[0] = dims[1] = dims[2] = parameters.mesh;

  assign_charge( system, &parameters, data, 0 );

  for(i=0;i<nm;i++)
    real_data[i] = data->Qmesh[2*i];

  write_mesh( out_file, real_data, dims, spacing, 1, "charge" );
}
