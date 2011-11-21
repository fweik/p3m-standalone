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
  FLOAT_TYPE last_error = 1e100;
  int j;
  int mesh_min, mesh_max, mesh_step=MESH_STEP, mesh_dir=1;
  int success_once = 0;
  int start_cao_run = 0;
  int cao_min, cao_max, cao_limit, cao_last;

  p_best.mesh = MESH_MAX;

  if( p->mesh != 0) {
    mesh_min = mesh_max = p->mesh;
  } else {
    mesh_min = MESH_MIN;
    mesh_max = MESH_MAX;
  }

  if( p->cao != 0 ) {
    cao_last = cao_min = cao_max = p->cao;
  } else {
    cao_min = CAO_MIN;
    cao_max = CAO_MAX;
    cao_last = CAO_MAX;
  }
  cao_start = CAO_MAX;

  it.rcut = p->rcut;

  TUNE_TRACE(printf("Starting tuning for '%s' with prec '%e'\n", m->method_name, precision);)
    TUNE_TRACE(printf("cao_min %d cao_max %d mesh_min %d mesh_max %d\n", cao_min, cao_max, mesh_min, mesh_max);)

    for(it.mesh = mesh_min; (it.mesh <= mesh_max) && (mesh_step >= 2); it.mesh+=mesh_dir*mesh_step ) {
      // If we were already successful with a smaller mesh
      // we can start at fastest cao so far.
    if( best_time < 1e200 ) 
      cao_start = (cao_last <= cao_min) ? cao_min : cao_last - 1;
    else
      cao_last = cao_start = cao_max;

    cao_limit = cao_min;

    if(start_cao_run == 0) {
      cao_start = cao_max;
      cao_limit = cao_max;
    }

    success = 0;
    printf("cao_start %d cao_limit %d\n", cao_start, cao_limit);
    for(it.cao = cao_start; (it.cao >= 2) && ( it.cao >= cao_limit ); it.cao--) {
      it.cao3 = it.cao * it.cao * it.cao;
      it.ip = it.cao - 1;

      Free_data(d);
      d = m->Init( s, &it );

      // reinit alpha loop
      direction = 1;
      alpha_step = ALPHA_STEP;
      last_error = min_error = 1e110;
      over_min = 0;

      for( it.alpha = 0.1; alpha_step >= 0.01; it.alpha += direction*alpha_step ) {

	error = m->Error( s, &it );
	//	TUNE_TRACE(printf("mesh %d cao %d rcut %e prec %e alpha %e dir %d overmin %d\n", it.mesh, it.cao, it.rcut, error, it.alpha, direction,over_min ););
	if( error <= precision && error < min_error ) {
	  min_error = error;
	  alpha_min_error = it.alpha;
	}

	if( (last_error - error) <= 0.0 ) {
	  direction = -direction;
          over_min = 1;
	}
	last_error=error;
	if( over_min == 1 )
	  alpha_step /= 2.0;
      }

       TUNE_TRACE(printf("fini mesh %d cao %d rcut %e prec %e alpha %e\n", it.mesh, it.cao, it.rcut, error, it.alpha ););
      
      if( min_error <= precision ) {
	it.alpha = alpha_min_error;   
	cao_last = it.cao;
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
    if((success == 1) && (start_cao_run == 0)) {
      //if we were successfull, reduce step and go back to see if 
      //accuracy can be matched with smaler mesh
      if( it.mesh > MESH_MIN ) {
	mesh_dir = -1;
	success_once = 1;
	mesh_step /= 2;
	if(mesh_step < 2)
	  mesh_step = 2;
      } else {
	start_cao_run = 1;
	mesh_dir = 1;
	mesh_step = MESH_STEP;
	cao_last = cao_max + 1;
      }
      
    } else {
      if( mesh_dir == -1) {
	start_cao_run = 1;
	mesh_dir = 1;
	mesh_step /= 2;
	cao_last = cao_max+1;
      if(mesh_step < 2)
	mesh_step = 2;
      if((it.mesh - mesh_step) <= mesh_step) {
	mesh_dir=1;
      }
      }else {
	mesh_dir = 1;
	if(start_cao_run == 1) {
	  mesh_step = 2;		       
	}
	else
	  mesh_step = MESH_STEP;
      }
    }
    if( ( success == 1 ) && ( time > 1.5*best_time ) ) {
      // No further gain expected with bigger mesh;
      break;
      } 
    if( it.mesh > 1.5*p_best.mesh )
      break;
  }
  if(d != NULL)
    Free_data(d);

  Free_forces(f);

<<<<<<< HEAD
=======
  if( success_once == 0 )
    return -1;

>>>>>>> b7944396ecc72170a59c39a616505b74c023fa2d
  *p = p_best;
  p->ip = p->cao - 1;
  p->cao3 = p->cao * p->cao * p->cao;

  TUNE_TRACE(printf("Using mesh %d cao %d rcut %e alpha %e with precision %e\n", p->mesh, p->cao, p->rcut, p->alpha, p->precision);)
    if( best_time < 1e100 )
      return 0;
    else
      return -1;
}  


