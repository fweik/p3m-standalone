#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "tuning.h"

#include "common.h"
#include "p3m-common.h"
#include "timings.h"
#include "interpol.h"
#include "realpart.h"

#define TUNE_DEBUG

#ifdef TUNE_DEBUG
  #include <stdio.h>
  #define TUNE_TRACE(A) A
#else
  #define TUNE_TRACE(A) 
#endif

parameters_t *Tune( const method_t *m, system_t *s, FLOAT_TYPE precision, FLOAT_TYPE rcut ) {
  parameters_t *p = Init_array( 1, sizeof( parameters_t ) );
  //Parameter iteraters, best parameter set
  parameters_t it, p_best;
  // Array to store forces
  forces_t *f = Init_forces(s->nparticles);

  data_t *d = NULL;

  it.rcut = rcut;

  FLOAT_TYPE best_time=1e250, time=1e240;

  FLOAT_TYPE last_error=1e100, error = -1.0;
  FLOAT_TYPE alpha_step = 0.1;
  int direction=1, over_min=0, success=0;

  TUNE_TRACE(printf("Starting tuning for '%s' with prec '%e'\n", m->method_name, precision);)

  for(it.mesh = MESH_MIN; it.mesh <= MESH_MAX; it.mesh+=MESH_STEP ) {
    TUNE_TRACE(printf("Trying mesh '%d'\n", it.mesh);)
    for(it.cao = CAO_MAX; it.cao >= CAO_MIN; it.cao--) {
      TUNE_TRACE(printf("Trying cao '%d'\n", it.cao);)
      it.cao3 = it.cao * it.cao * it.cao;
      it.ip = it.cao - 1;

      Free_data(d);
      d = m->Init( s, &it );

      // reinit alpha loop
      direction = 1;
      alpha_step = ALPHA_STEP;

      for( it.alpha = ALPHA_STEP_MIN; alpha_step <= ALPHA_STEP_MIN; it.alpha += direction*alpha_step ) {
	TUNE_TRACE(printf("alpha %e\n", it.alpha);)
	error = m->Error( s, &it );
	if( (last_error - error) >= 1e-3 ) {
	  continue;
	}
	if( (last_error - error) <= 1e-3 ) {
	  direction = -direction;
          over_min = 1;
	}
	if( over_min == 1 )
	  alpha_step /= 2.0;
      } 
      /*
       TUNE_TRACE(puts("Starting timing...");)

	it.rcut = last_success;   
	success = 1;

	Init_neighborlist( s, &it, d );

	m->Influence_function( s, &it, d );

	time = MPI_Wtime();
	Calculate_forces ( m, s, &it, d, f );
	time = MPI_Wtime() - time;

	Free_neighborlist(d);

	TUNE_TRACE(printf("\n mesh %d cao %d rcut %e time %e prec %e alpha %e\n", it.mesh, it.cao, it.rcut, time, error, it.alpha );)

	if( time < best_time ) {
	  p_best = it;
	  p_best.precision = error;
	  best_time = time;
	  TUNE_TRACE(printf("New best. mesh %d cao %d rcut %e time %e prec %e alpha %e\n", p_best.mesh, p_best.cao, p_best.rcut, time, error, p_best.alpha );)         
	}
      } else {
	// If we didn't get the precision here, there is no point going to lower cao.
	break;
	} */
    }
    if( time > 1.2*best_time ) {
      // No further gain expected with bigger mesh;
      break;
    }
  }
  if(d != NULL)
    Free_data(d);

  Free_forces(f);

  if( success == 0 )
    return NULL;

  *p = p_best;
  p->ip = p->cao - 1;
  p->cao3 = p->cao * p->cao * p->cao;

  TUNE_TRACE(printf("Using mesh %d cao %d rcut %e alpha %e with precision %e\n", p->mesh, p->cao, p->rcut, p->alpha, p->precision);)

  return p;
}  


