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

int Tune( const method_t *m, system_t *s, parameters_t *p, FLOAT_TYPE precision ) {
  //Parameter iteraters, best parameter set
  parameters_t it, p_best;
  // Array to store forces
  forces_t *f = Init_forces(s->nparticles);

  data_t *d = NULL;

  FLOAT_TYPE best_time=1e250, time=1e240;

  FLOAT_TYPE min_error=1e110, error = -1.0;
  FLOAT_TYPE alpha_step = 0.1;
  FLOAT_TYPE alpha_min_error = 0.0;
  int direction=1, over_min=0, success=0, cao_start;
  FLOAT_TYPE last_error = 0.0;

  int mesh_min, mesh_max;
  int cao_min, cao_max;

  if( p->mesh != 0) {
    mesh_min = mesh_max = p->mesh;
  } else {
    mesh_min = MESH_MIN;
    mesh_max = MESH_MAX;
  }

  if( p->cao != 0 ) {
    cao_min = cao_max = p->cao;
  } else {
    cao_min = CAO_MIN;
    cao_max = CAO_MAX;
  }
  cao_start = CAO_MAX;

  TUNE_TRACE(printf("Starting tuning for '%s' with prec '%e'\n", m->method_name, precision);)
    TUNE_TRACE(printf("cao_min %d cao_max %d mesh_min %d mesh_max %d\n", cao_min, cao_max, mesh_min, mesh_max);)

  for(it.mesh = mesh_min; it.mesh <= mesh_max; it.mesh+=MESH_STEP ) {
      // If we were already successful with a smaller mesh
      // we can start at fastest cao so far.
    if( best_time < 1e200 ) 
      cao_start = (p_best.cao <= cao_min) ? cao_min : p_best.cao - 1;
    else
      cao_start = cao_max;

    for(it.cao = cao_start; (it.cao >= 2) && ( it.cao >= cao_min ); it.cao--) {

      it.cao3 = it.cao * it.cao * it.cao;
      it.ip = it.cao - 1;

      Free_data(d);
      d = m->Init( s, &it );

      // reinit alpha loop
      direction = 1;
      alpha_step = ALPHA_STEP;
      min_error = 1e110;

      for( it.alpha = 0.2; alpha_step >= 0.01; it.alpha += direction*alpha_step ) {
	if( p->rcut != 0.0 ) {
	  it.rcut = p->rcut;
	} else {

	}

	error = m->Error( s, &it );

	if( error <= precision && error < min_error ) {
	  min_error = error;
	  alpha_min_error = it.alpha;
	}
	if( (last_error - error) >= 0.0 ) {
	  last_error = error;
	  continue;
	}
	if( (last_error - error) <= 0.0 ) {
	  last_error=error;
	  direction = -direction;
          over_min = 1;
	}
	if( over_min == 1 )
	  alpha_step /= 2.0;
      }

      TUNE_TRACE(printf("mesh %d cao %d rcut %e prec %e alpha %e\n", it.mesh, it.cao, it.rcut, error, it.alpha ););
      
      if( min_error <= precision ) {
	it.alpha = alpha_min_error;   
	success = 1;
      } else {
	break;
      }

      if( min_error > precision )
	break;

      TUNE_TRACE(puts("Starting timing..."););

      m->Influence_function( s, &it, d );

      time = MPI_Wtime();
      m->Kspace_force( s, &it, d, f );
      time = MPI_Wtime() - time;

      TUNE_TRACE(printf("\n mesh %d cao %d rcut %e time %e prec %e alpha %e\n", it.mesh, it.cao, it.rcut, time, error, it.alpha ););
	
      if( time < best_time ) {
	p_best = it;
	p_best.precision = error;
	best_time = time;
	TUNE_TRACE(printf("New best. mesh %d cao %d rcut %e time %e prec %e alpha %e\n", p_best.mesh, p_best.cao, p_best.rcut, time, error, p_best.alpha );)         ;
      }
      else {
	// If we didn't get the precision here, there is no point going to lower cao.
	break;
      } 
    }
    if( ( time > 1.1*best_time ) || ( it.mesh > 2*p_best.mesh ) ) {
      // No further gain expected with bigger mesh;
      break;
    }
  }
  if(d != NULL)
    Free_data(d);

  Free_forces(f);

  if( success == 0 )
    return -1;

  *p = p_best;
  p->ip = p->cao - 1;
  p->cao3 = p->cao * p->cao * p->cao;

  TUNE_TRACE(printf("Using mesh %d cao %d rcut %e alpha %e with precision %e\n", p->mesh, p->cao, p->rcut, p->alpha, p->precision);)
    return 0;
}  


