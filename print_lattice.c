#include "types.h"
#include "io.h"
#include "interpol.h"
#include "charge-assign.h"
#include "p3m-common.h"
#include "common.h"

int main(int argc, char **argv) {

  cmd_parameters_t params = { NULL, 0, NULL, 0 };
  char *out_file, *pos_file;

  parameters_t parameters;
  system_t *system;
  data_t *data;
  method_t m = { 255, "assign only", METHOD_FLAG_Qmesh | METHOD_FLAG_ca, NULL, NULL, NULL, NULL };

  double *real_data;
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

  dims[0] = dims[1] = dims[2] = parameters.mesh;

  assign_charge( system, &parameters, data, 0 );

  for(i=0;i<nm;i++)
    real_data[i] = data->Qmesh[2*i];

  write_mesh( out_file, real_data, dims, spacing, 1, "charge" );
}
