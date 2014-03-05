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
#include <float.h>

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

const int smooth_numbers[] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 36, 40, 42, 44, 48, 50, 52, 54, 56, 60, 64, 66, 70, 72, 78, 80, 84, 88, 90, 96, 98, 100, 104, 108, 110, 112, 120, 126, 128, 130, 132, 140, 144, 150, 154, 156, 160, 162, 168, 176, 180, 182, 192, 196, 198, 200, 208, 210, 216, 220, 224, 234, 240, 242, 250, 252, 256, 260, 264, 270, 280, 288, 294, 300};

const int smooth_numbers_n = sizeof(smooth_numbers)/sizeof(int);

typedef struct {
  double *mesh_timings;
  runtime_t *cao_timings;
  int n;
} parameter_timings_t;

parameter_timings_t *pt = NULL;
int n_pt = 0;

static int
compare_ints (const void *a, const void *b)
{
  return (*(int *)(a) - *(int *)(b));
}

double time_mesh(const method_t  *m, system_t *s, parameters_t *p) {
parameters_t mp = *p;
mp.cao = CAO_MIN;
mp.cao3 = mp.cao*mp.cao*mp.cao;
mp.ip = mp.cao - 1;
mp.tuning = 1;
data_t *d = m->Init(s, &mp);
double res  = 0.0;

  m->Kspace_force( s, &mp, d, s->reference );

  res = d->runtime.t_g;

  Free_data(d);

  return res;
}

runtime_t time_full(const method_t  *m, system_t *s, parameters_t *p) {
parameters_t mp = *p;
mp.tuning = 1;
data_t *d = m->Init(s, &mp);
runtime_t ret;

  m->Kspace_force( s, &mp, d, s->reference );

  ret = d->runtime;

  Free_data(d);

  return ret;
}


runtime_t time_cao(const method_t *m, system_t *s, parameters_t *p) {
  parameters_t mp = *p;
  mp.tuning = 1;
  data_t *d = m->Init(s, &mp);
  runtime_t res;

  m->Kspace_force( s, &mp, d, s->reference );

  res = d->runtime;
  res.t_g = 0.0;

/* TUNE_TRACE(printf("time_cao(%d): t_c = %e, t_f = %e\n", p->cao, d->runtime.t_c, d->runtime.t_f)); */

  Free_data(d);

  return res;
}

runtime_t get_timing(const method_t *m, system_t *s, parameters_t *p) {
  /* printf("get_timing(method_id %d, mesh %d, cao %d pt %p)\n", m->method_id, p->mesh, p->cao, pt); */

runtime_t ret;

  /* if(n_pt == 0) { */
  /*   pt = Init_array(1, sizeof(parameter_timings_t)); */
  /*   n_pt = 1; */
  /*   pt[0].mesh_timings = NULL; */
  /*   pt[0].cao_timings = NULL; */
  /* } */
  /* if(m->method_id >= n_pt) { */
  /*   pt = realloc(pt, (m->method_id+1) * sizeof(parameter_timings_t)); */

  /*   for(int i = n_pt-1; i <= m->method_id; i++) { */
  /*     pt[i].mesh_timings = NULL; */
  /*     pt[i].cao_timings = NULL; */
  /*   } */
  /*   n_pt = m->method_id+1; */
  /* } */

  /* if(pt[m->method_id].mesh_timings == NULL) { */
  /*   pt[m->method_id].mesh_timings = Init_array(smooth_numbers_n, sizeof(double)); */
  /*   memset(pt[m->method_id].mesh_timings, 0, smooth_numbers_n * sizeof(double)); */
  /* } */
  /* if((pt[m->method_id].cao_timings == NULL) || (pt[m->method_id].n != s->nparticles)) { */
  /*   puts("Resetting cao timings."); */
  /*   Free_array(pt[m->method_id].cao_timings); */
  /*   pt[m->method_id].cao_timings = Init_array(CAO_MAX+1, sizeof(runtime_t)); */
  /*   memset(pt[m->method_id].cao_timings, 0, (CAO_MAX+1) * sizeof(runtime_t)); */
  /*   pt[m->method_id].n = s->nparticles; */
  /* }   */

  /* int *mesh_of = bsearch(&(p->mesh), smooth_numbers, smooth_numbers_n, sizeof(int), compare_ints); */
  /* int mesh_id = mesh_of - smooth_numbers; */

  /* if(pt[m->method_id].mesh_timings[mesh_id] <= 0) { */
  /*   pt[m->method_id].mesh_timings[mesh_id] = time_mesh(m,s,p); */
  /*   /\* TUNE_TRACE(printf("mesh miss %d (id %d), time %e\n", p->mesh, mesh_id, pt[m->method_id].mesh_timings[mesh_id]);); *\/ */
  /* } */
  /* if(pt[m->method_id].cao_timings[p->cao].t_c <= 0) { */
  /*   pt[m->method_id].cao_timings[p->cao] = time_cao(m,s,p); */
  /*   /\* TUNE_TRACE(printf("cao miss %d, time %e\n", p->cao, s->nparticles * pt[m->method_id].cao_timings[p->cao]);); *\/ */
  /* } */

ret = time_full(m, s, p);

//  ret = pt[m->method_id].cao_timings[p->cao];
//  ret.t_g = pt[m->method_id].mesh_timings[mesh_id];


  ret.t = ret.t_c + ret.t_g + ret.t_f;

  return  ret;
}

runtime_t Tune( const method_t *m, system_t *s, parameters_t *p, FLOAT_TYPE precision ) {
  //Parameter iteraters, best parameter set
  parameters_t it, p_best;
  // Array to store forces
// @TODO: Make copy of reference forces.

  forces_t *f = Init_forces(s->nparticles);

  it.prefactor = 1.0;

  runtime_t best_time;
  best_time.t = DBL_MAX;
  runtime_t time;

  FLOAT_TYPE error = -1.0;
  FLOAT_TYPE V = pow( s->length, 3);
  int cao_start;
  int mesh_it_min=0, mesh_it_max=smooth_numbers_n, mesh_it;
  int success_once = 0;
  int cao_min, cao_max, cao_limit, cao_last;

  runtime_t ret;
  ret.t = -1;

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
  it.tuning = p->tuning;

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

    /* TUNE_TRACE(printf("cao_start %d cao_limit %d\n", cao_start, cao_limit);); */

    it.cao = cao_min;
    it.cao3 = it.cao * it.cao * it.cao;
    it.ip = it.cao - 1;

    // Check if we allready are slower than the best timing.
    // Then there is no point in going on, it will only get worse.
    runtime_t min_time;
    min_time = get_timing(m, s, &it);

    if( min_time.t > best_time.t) {
      TUNE_TRACE(printf("Best possible time for (%d %d) = %e slower than best %e\n", it.mesh, it.cao, min_time.t, best_time.t););
      break;
    }


    for(it.cao = cao_start; (it.cao >= cao_min) && ( it.cao >= cao_limit ); it.cao--) {
      it.cao3 = it.cao * it.cao * it.cao;
      it.ip = it.cao - 1;

      error = m->Error( s, &it);

      /* TUNE_TRACE(printf("fini mesh %d cao %d rcut %e prec %e alpha %e\n", it.mesh, it.cao, it.rcut, error, it.alpha );); */
      
      if( error <= precision ) {
	cao_last = it.cao;
	success_once = 1;
      } else {
	break;
      }

      if( (it.mesh >= p_best.mesh) && (it.cao >= p_best.cao))
	break;

      /* TUNE_TRACE(puts("Starting timing...");); */

      time = get_timing(m, s, &it);
     
      TUNE_TRACE(printf("\n Timing mesh %d cao %d rcut %e time %e prec %e alpha %e\n", it.mesh, it.cao, it.rcut, time.t, error, it.alpha ););
	
      if( time.t < best_time.t ) {
	p_best = it;
	p_best.precision = error;
	best_time = time;
	ret = best_time;
	TUNE_TRACE(printf("New best. mesh %d cao %d rcut %e time %e prec %e alpha %e\n", p_best.mesh, p_best.cao, p_best.rcut, time.t, error, p_best.alpha );)         ;
      }
    }
  }
  Free_forces(f);

  if( success_once == 0 )
    return ret;

  *p = p_best;
  p->ip = p->cao - 1;
  p->cao3 = p->cao * p->cao * p->cao;

  TUNE_TRACE(printf("Using mesh %d cao %d rcut %e alpha %e with precision %e time %e\n", p->mesh, p->cao, p->rcut, p->alpha, p->precision, best_time.t);)
  return ret;
}


