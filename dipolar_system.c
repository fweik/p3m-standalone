#include "types.h"
#include "io.h"
#include "generate_system.h"
#include "common.h"
#include "error.h"
#include "p3m-ad-i.h"
#include "p3m-ik.h"
#include "p3m-ad.h"
#include "p3m-ad-self-forces.h"

#include <math.h>

int main(void) {
  system_t *system = generate_system( SYSTEM_RANDOM, 200, 10.0, 1.0);
  parameters_t parameters;
  forces_t *forces;
  FLOAT_TYPE ref_prec = -1.0;
  method_t method = method_p3m_ad;
  data_t *data;
  error_t error;
  FILE *f = fopen("out.dat", "w");

  parameters.prefactor = 1.0;
  parameters.rcut = 3.0;
  parameters.cao = 7;
  parameters.ip = parameters.cao - 1;
  parameters.cao3 =  parameters.cao* parameters.cao* parameters.cao;
  parameters.mesh = 16;

  forces = Init_forces(system->nparticles);
 
  ref_prec = Calculate_reference_forces( system, &parameters );

  printf("Reference precision (Ewald) %e\n", ref_prec);

  data = method.Init(system, &parameters);

  for ( parameters.alpha=1.0; parameters.alpha<=1.0; parameters.alpha+=0.1 ) {
    method.Influence_function ( system, &parameters, data );

    Init_self_forces( system, &parameters, data );

    Calculate_forces ( &method, system, &parameters, data, forces );

    Substract_self_forces(system, &parameters, data, forces );

    error = Calculate_errors ( system, forces );

    if ( method.Error != NULL ) {
      double estimate =  method.Error ( system, &parameters );
      printf ( "%8lf\t%8e\t%8e\t %8e %8e\n", parameters.alpha, error.f / sqrt(system->nparticles) , estimate,
	       error.f_r / sqrt(system->nparticles), error.f_k / sqrt(system->nparticles));
      fprintf ( f, "%8lf\t%8e\t%8e\n", parameters.alpha, error.f / sqrt(system->nparticles), estimate);
      fflush(stdout);
      fflush(f);
    }
  }


  return 0;
}
