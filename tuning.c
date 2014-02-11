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

#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#include "tuning.h"

#include "common.h"
#include "p3m-common.h"
#include "timings.h"
#include "interpol.h"
#include "realpart.h"

static int dummy_data_initialized = 0;
static data_t dummy_ad_complex;
static data_t dummy_ad_real;

static data_t dummy_ik_complex;
static data_t dummy_ik_real;


const int smooth_numbers[] = {4, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 36, 40, 42, 44, 48, 50, 52, 54, 56, 60, 64, 66, 70, 72, 78, 80, 84, 88, 90, 96, 98, 100, 104, 108, 110, 112, 120, 126, 128, 130, 132, 140, 144, 150, 154, 156, 160, 162, 168, 176, 180, 182, 192, 196, 198, 200, 208, 210, 216, 220, 224, 234, 240, 242, 250, 252, 256, 260, 264, 270, 280, 288, 294, 300};

const int smooth_numbers_n = sizeof(smooth_numbers)/sizeof(int);

void Init_dummy(int max_part) {
  int mesh = smooth_numbers[smooth_numbers_n-1];
  int mesh3 = mesh*mesh*mesh ;

  data_t *d = Init_array(1, sizeof(data_t));

  d->mesh = mesh;

  d->Qmesh = Init_array(2*mesh3, sizeof(FLOAT_TYPE));

  d->Fmesh = Init_vector_array(2*mesh3);
  d->Dn = Init_array(d->mesh, sizeof(FLOAT_TYPE));
  Init_differential_operator(d);

  d->nshift = Init_array(d->mesh, sizeof(FLOAT_TYPE));
  Init_nshift(d);


  int max = 2;

  for (int i = 0; i < max; i++) {
    d->dQdx[i] = Init_array( max_part*7*7*7, sizeof(FLOAT_TYPE) );
    d->dQdy[i] = Init_array( max_part*7*7*7, sizeof(FLOAT_TYPE) );
    d->dQdz[i] = Init_array( max_part*7*7*7, sizeof(FLOAT_TYPE) );
  }
      
  for (int i = 0; i < max; i++) {
    d->cf[i] = Init_array( 7*7*7 * max_part, sizeof(FLOAT_TYPE));
    d->ca_ind[i] = Init_array( 3*max_part, sizeof(int));
  }
	
  d->inter = Init_interpolation( 6, 1 );
	
  d->G_hat = Init_array(mesh3, sizeof(FLOAT_TYPE));

  d->forward_plans = 1;
  d->backward_plans = 3;

  memcpy(&dummy_ad_real, d, sizeof(data_t));
  memcpy(&dummy_ad_complex, d, sizeof(data_t));
  memcpy(&dummy_ik_real, d, sizeof(data_t));
  memcpy(&dummy_ik_complex, d, sizeof(data_t));

  int l;

  dummy_ik_complex.forward_plan[0] = FFTW_PLAN_DFT_3D ( mesh, mesh, mesh, ( FFTW_COMPLEX * ) dummy_ik_complex.Qmesh, ( FFTW_COMPLEX * ) dummy_ik_complex.Qmesh, FFTW_FORWARD, FFTW_PATIENT );

  for ( l=0;l<3;l++ ) {
    dummy_ik_complex.backward_plan[l] = FFTW_PLAN_DFT_3D ( mesh, mesh, mesh, ( FFTW_COMPLEX * ) ( dummy_ik_complex.Fmesh->fields[l] ), ( FFTW_COMPLEX * ) ( dummy_ik_complex.Fmesh->fields[l] ), FFTW_BACKWARD, FFTW_PATIENT );
  }

    dummy_ad_complex.forward_plan[0] = FFTW_PLAN_DFT_3D ( mesh, mesh, mesh, ( FFTW_COMPLEX * ) dummy_ad_complex.Qmesh, ( FFTW_COMPLEX * ) dummy_ad_complex.Qmesh, FFTW_FORWARD, FFTW_PATIENT );

    dummy_ad_complex.backward_plan[0] = FFTW_PLAN_DFT_3D ( mesh, mesh, mesh, ( FFTW_COMPLEX * ) ( dummy_ad_complex.Qmesh ), ( FFTW_COMPLEX * ) ( dummy_ad_complex.Qmesh ), FFTW_BACKWARD, FFTW_PATIENT );

    dummy_ik_real.forward_plan[0] = FFTW_PLAN_DFT_R2C_3D ( mesh, mesh, mesh, dummy_ik_real.Qmesh, (FFTW_COMPLEX *)dummy_ik_real.Qmesh, FFTW_PATIENT );

    for ( l=0;l<3;l++ ) {
        dummy_ik_real.backward_plan[l] = FFTW_PLAN_DFT_C2R_3D ( mesh, mesh, mesh, ( FFTW_COMPLEX * ) ( dummy_ik_real.Fmesh->fields[l] ), ( dummy_ik_real.Fmesh->fields[l] ), FFTW_PATIENT );
    }

    dummy_ad_real.forward_plan[0] = FFTW_PLAN_DFT_R2C_3D ( mesh, mesh, mesh, dummy_ad_real.Qmesh, ( FFTW_COMPLEX * ) dummy_ad_real.Qmesh, FFTW_PATIENT );

    dummy_ad_real.backward_plan[0] = FFTW_PLAN_DFT_C2R_3D ( mesh, mesh, mesh, ( FFTW_COMPLEX * ) ( dummy_ad_real.Qmesh ),  dummy_ad_real.Qmesh,  FFTW_PATIENT );

  dummy_data_initialized = 1;
}


//#define TUNE_DEBUG

#ifdef TUNE_DEBUG
  #include <stdio.h>
  #define TUNE_TRACE(A) A
#else
  #define TUNE_TRACE(A) 
#endif

timing_t Tune( const method_t *m, system_t *s, parameters_t *p, FLOAT_TYPE precision ) {
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

  timing_t ret = { -1, 0, N_TUNING_SAMPLES };

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

      if(!dummy_data_initialized)
d = m->Init( s, &it );
      else if ( m->method_id == 0) {
	d = &dummy_ik_complex;
      } else {
	d = m->Init( s, &it );
      }
      
      
      TUNE_TRACE(puts("Starting timing..."););

      double avg = 0, sgm = 0;
      for(int i = 0; i < N_TUNING_SAMPLES; i++) {
      time = MPI_Wtime();

      m->Kspace_force( s, &it, d, f );

      time = MPI_Wtime() - time;
      avg += time;
      sgm += time*time;
      }

      avg /= N_TUNING_SAMPLES;
      sgm = sqrt(sgm/N_TUNING_SAMPLES - avg*avg);

      TUNE_TRACE(printf("\n mesh %d cao %d rcut %e time %e prec %e alpha %e\n", it.mesh, it.cao, it.rcut, time, error, it.alpha ););
	
      if( time < best_time ) {
	p_best = it;
	p_best.precision = error;
	best_time = time;
	ret.avg = avg;
	ret.sgm = sgm;
	ret.n = N_TUNING_SAMPLES;
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
    return ret;

  *p = p_best;
  p->ip = p->cao - 1;
  p->cao3 = p->cao * p->cao * p->cao;

  TUNE_TRACE(printf("Using mesh %d cao %d rcut %e alpha %e with precision %e time %e\n", p->mesh, p->cao, p->rcut, p->alpha, p->precision, best_time);)
  return ret;
}


