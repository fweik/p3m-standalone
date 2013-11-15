#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "tuning.h"

#include "common.h"
#include "p3m-common.h"
#include "timings.h"
#include "interpol.h"
#include "realpart.h"

const int smooth_numbers[] = {4, 5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 36, 40, 42, 44, 48, 50, 52, 54, 56, 60, 64, 66, 70, 72, 78, 80, 84, 88, 90, 96, 98, 100, 104, 108, 110, 112, 120, 126, 128, 130, 132, 140, 144, 150, 154, 156, 160, 162, 168, 176, 180, 182, 192, 196, 198, 200, 208, 210, 216, 220, 224, 234, 240, 242, 250, 252, 256, 260, 264, 270, 280, 288, 294, 300};

const int smooth_numbers_n = 122;

//#define TUNE_DEBUG

#ifdef TUNE_DEBUG
  #include <stdio.h>
  #define TUNE_TRACE(A) A
#else
  #define TUNE_TRACE(A) 
#endif

FLOAT_TYPE Tune( const method_t *m, system_t *s, parameters_t *p, FLOAT_TYPE precision ) {
  //Parameter iteraters, best parameter set
  parameters_t it, p_best;
  // Array to store forces
  forces_t *f = Init_forces(s->nparticles);

  data_t *d = NULL;

  it.prefactor = 1.0;

  FLOAT_TYPE best_time=1e250, time=1e240;

  FLOAT_TYPE error = -1.0;
  FLOAT_TYPE V = pow( s->length, 3);
  int cao_start;
  int mesh_it_min=0, mesh_it_max=smooth_numbers_n, mesh_it;
  int success_once = 0;
  int cao_min, cao_max, cao_limit, cao_last;

  p_best.mesh = smooth_numbers[mesh_it_max] + 1;
  p_best.cao = CAO_MAX + 1;

  if( p->cao != 0 ) {
    cao_last = cao_min = cao_max = p->cao;
  } else {
    cao_min = CAO_MIN;
    cao_max = CAO_MAX;
    cao_last = CAO_MAX;
  }
  cao_start = CAO_MAX;

  it.rcut = p->rcut;

  TUNE_TRACE(printf("Starting tuning for '%s' with prec '%e'\n", m->method_name, precision););
  TUNE_TRACE(printf("cao_min %d cao_max %d\n", cao_min, cao_max););

  it.alpha = SQRT(-LOG((precision*SQRT(s->nparticles*it.rcut*V))/(2*SQRT(2)*s->q2)))/it.rcut;

  for(mesh_it = mesh_it_min; mesh_it <= mesh_it_max; mesh_it++ ) {
    it.mesh = smooth_numbers[mesh_it];
    if( success_once == 1 ) {
      if(cao_last <= cao_min) {
	break;
      } else {
	cao_start = cao_last - 1;
      }
    } else {
      cao_last = cao_start = cao_max;
    }

    cao_limit = cao_min;

    TUNE_TRACE(printf("cao_start %d cao_limit %d\n", cao_start, cao_limit););

    for(it.cao = cao_start; (it.cao >= cao_min) && ( it.cao >= cao_limit ); it.cao--) {
      it.cao3 = it.cao * it.cao * it.cao;
      it.ip = it.cao - 1;

      error = m->Error( s, &it);

      TUNE_TRACE(printf("fini mesh %d cao %d rcut %e prec %e alpha %e\n", it.mesh, it.cao, it.rcut, error, it.alpha ););
      
      if( error <= precision ) {
	cao_last = it.cao;
	success_once = 1;
      } else {
	break;
      }

      if( (it.mesh >= p_best.mesh) && (it.cao >= p_best.cao))
	break;

      TUNE_TRACE(puts("Initializing..."););

      Free_data(d);
      d = m->Init( s, &it );

      TUNE_TRACE(puts("Starting timing..."););

      const int samples = 100;
      double avg = 0;

      for(int i = 0; i < samples; i++ ) {
	time = MPI_Wtime();
	m->Kspace_force( s, &it, d, f );
        time = MPI_Wtime() - time;
	/* printf("run %d time %lf\n", i+1, time);  */
        avg += time;
      }
      avg /= samples;
      time = avg;

      TUNE_TRACE(printf("\n mesh %d cao %d rcut %e time %e prec %e alpha %e\n", it.mesh, it.cao, it.rcut, time, error, it.alpha ););
	
      if( time < best_time ) {
	p_best = it;
	p_best.precision = error;
	best_time = time;
	TUNE_TRACE(printf("New best. mesh %d cao %d rcut %e time %e prec %e alpha %e\n", p_best.mesh, p_best.cao, p_best.rcut, time, error, p_best.alpha );)         ;
      }
    }
    if( (success_once == 1) && (it.mesh > (p_best.mesh + 10)) )
      break;
  }
  if(d != NULL)
    Free_data(d);

  Free_forces(f);

  if( success_once == 0 )
    return -1;

  *p = p_best;
  p->ip = p->cao - 1;
  p->cao3 = p->cao * p->cao * p->cao;

  TUNE_TRACE(printf("Using mesh %d cao %d rcut %e alpha %e with precision %e time %e\n", p->mesh, p->cao, p->rcut, p->alpha, p->precision, best_time);)
    if( best_time < 1e100 )
      return best_time;
    else
      return -1;
}


