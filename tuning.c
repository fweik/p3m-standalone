/**    Copyright (C) 2011,2012,2013,2014 Florian Weik <fweik@icp.uni-stuttgart.de>

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
#include <string.h>
#include <float.h>

#include "tuning.h"

#include "common.h"
#include "p3m-common.h"
#include "interpol.h"
#include "realpart.h"

//#define TUNE_DEBUG

#ifdef TUNE_DEBUG
  #include <stdio.h>
  #define TUNE_TRACE(A) A
#else
  #define TUNE_TRACE(A) 
#endif

const int smooth_numbers[] = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 36, 40, 42, 44, 48, 50, 52, 54, 56, 60, 64, 66, 70, 72, 78, 80, 84, 88, 90, 96, 98, 100, 104, 108, 110, 112, 120, 126, 128, 130, 132, 140, 144, 150, 154, 156, 160, 162, 168, 176, 180, 182, 192, 196, 198, 200, 208, 210, 216, 220, 224, 234, 240, 242, 250, 252, 256, 260, 264, 270, 280, 288, 294, 300};
const int smooth_numbers_n = sizeof(smooth_numbers)/sizeof(int);

const int powers_of_two[] = { 2, 2 << 1, 2 << 2,2 << 3,2 << 4,2 << 5,2 << 6,2 << 7,2 << 8,2 << 9,2 << 10};
const int powers_of_two_n = sizeof(powers_of_two)/sizeof(int);

double time_series[N_TUNING_SAMPLES];

int time_fluctuation_hist[1000];
int n_hist = 0;


typedef struct {
  timing_t *mesh_timings;
  runtime_stat_t *cao_timings;
  int n;
} parameter_timings_t;

parameter_timings_t *pt = NULL;
int n_pt = 0;

void time_hist(double v) {
  const int bins = sizeof(time_fluctuation_hist) / sizeof(int);
  double x = fabs(v);
  int bin = (int)floor(x*bins);
  if (bin >= bins)
    bin = bins-1;

  time_fluctuation_hist[bin]++;
  n_hist++;
}

timing_t time_mesh(const method_t  *m, system_t *s, parameters_t *p) {
parameters_t mp = *p;
mp.cao = CAO_MIN;
mp.cao3 = mp.cao*mp.cao*mp.cao;
mp.ip = mp.cao - 1;
mp.tuning = 1;
data_t *d = m->Init(s, &mp);
timing_t res  = {0.0, 0.0, 1e99, 0.0, N_TUNING_SAMPLES};

 for(int i = 0; i < N_TUNING_SAMPLES; i++) {
   m->Kspace_force( s, &mp, d, s->reference );

   res.avg += d->runtime.t_g;
   res.min = ( d->runtime.t_g < res.min ) ? d->runtime.t_g : res.min;
   res.max = ( d->runtime.t_g > res.max ) ? d->runtime.t_g : res.max;
 }

 res.avg /= res.n;

 Free_data(d);
 return res;
}

runtime_stat_t time_full(const method_t  *m, system_t *s, parameters_t *p) {
  parameters_t mp = *p;
  mp.tuning = 1;
  data_t *d = m->Init(s, &mp);
  runtime_stat_t res  = {{0.0, 0.0, 1e99, 0.0, N_TUNING_SAMPLES},
			 {0.0, 0.0, 1e99, 0.0, N_TUNING_SAMPLES},
			 {0.0, 0.0, 1e99, 0.0, N_TUNING_SAMPLES},
			 {0.0, 0.0, 1e99, 0.0, N_TUNING_SAMPLES}};


  for(int i = 0; i < N_TUNING_SAMPLES; i++) {

    m->Kspace_force( s, &mp, d, s->reference );

    res.t_g.avg += d->runtime.t_g;
    res.t_g.min = ( d->runtime.t_g < res.t_g.min ) ? d->runtime.t_g : res.t_g.min;
    res.t_g.max = ( d->runtime.t_g > res.t_g.max ) ? d->runtime.t_g : res.t_g.max;

    res.t_f.avg += d->runtime.t_f;
    res.t_f.min = ( d->runtime.t_f < res.t_f.min ) ? d->runtime.t_f : res.t_f.min;
    res.t_f.max = ( d->runtime.t_f > res.t_f.max ) ? d->runtime.t_f : res.t_f.max;

    res.t_c.avg += d->runtime.t_c;
    res.t_c.min = ( d->runtime.t_c < res.t_c.min ) ? d->runtime.t_c : res.t_c.min;
    res.t_c.max = ( d->runtime.t_c > res.t_c.max ) ? d->runtime.t_c : res.t_c.max;

    time_series[i] = d->runtime.t_c + d->runtime.t_g + d->runtime.t_f;

    memset(&(d->runtime), 0, sizeof(runtime_t));
  }

  res.t_f.avg /= res.t_f.n;
  res.t_c.avg /= res.t_c.n;
  res.t_g.avg /= res.t_g.n;

  res.t.avg = res.t_c.avg + res.t_g.avg + res.t_f.avg;
  res.t.min = res.t_c.min + res.t_g.min + res.t_f.min;
  res.t.max = res.t_c.max + res.t_g.max + res.t_f.max;

  if(0) {
    double ep = (res.t.max - res.t.avg) / res.t.avg;
    double em = (res.t.avg - res.t.min) / res.t.avg;

    printf("warning rel timing fluctuations +%e -%e\n", ep, em);
    printf("warning time series (avg %e)\n", res.t.avg);
    for(int i = 0; i < N_TUNING_SAMPLES; i++) {
      printf("warning %d %e %e\n", i, time_series[i], (time_series[i] - res.t.avg)/res.t.avg);
      time_hist((time_series[i] - res.t.avg)/res.t.avg);
    }
    puts("warning");
  }
  
  Free_data(d);

  return res;
}


void write_hist(void) {
  puts("write_host()");
  FILE *f = fopen("hist.dat", "w");
  const int bins = sizeof(time_fluctuation_hist) / sizeof(int);

  for(int i = 0; i < bins; i++) {
    fprintf(f, "%d %e\n", i, (double)(time_fluctuation_hist[i]));
  }
  fclose(f);
}

runtime_stat_t time_cao(const method_t *m, system_t *s, parameters_t *p) {
  parameters_t mp = *p;
  mp.tuning = 1;
  data_t *d = m->Init(s, &mp);
  runtime_stat_t res  = {{0.0, 0.0, 1e99, 0.0, N_TUNING_SAMPLES},
			 {0.0, 0.0, 1e99, 0.0, N_TUNING_SAMPLES},
			 {0.0, 0.0, 1e99, 0.0, N_TUNING_SAMPLES},
			 {0.0, 0.0, 1e99, 0.0, N_TUNING_SAMPLES}};

 for(int i = 0; i < N_TUNING_SAMPLES; i++) {

   m->Kspace_force( s, &mp, d, s->reference );

   res.t_f.avg += d->runtime.t_f;
   res.t_f.min = ( d->runtime.t_f < res.t_f.min ) ? d->runtime.t_f : res.t_f.min;
   res.t_f.max = ( d->runtime.t_f > res.t_f.max ) ? d->runtime.t_f : res.t_f.max;

   res.t_c.avg += d->runtime.t_c;
   res.t_c.min = ( d->runtime.t_c < res.t_c.min ) ? d->runtime.t_c : res.t_c.min;
   res.t_c.max = ( d->runtime.t_c > res.t_c.max ) ? d->runtime.t_c : res.t_c.max;
 }

 res.t_f.avg /= res.t_f.n;
 res.t_c.avg /= res.t_c.n;
 
/* TUNE_TRACE(printf("time_cao(%d): t_c = %e, t_f = %e\n", p->cao, d->runtime.t_c, d->runtime.t_f)); */

  Free_data(d);

  return res;
}

runtime_stat_t get_timing(const method_t *m, system_t *s, parameters_t *p) {
  /* printf("get_timing(method_id %d, mesh %d, cao %d pt %p)\n", m->method_id, p->mesh, p->cao, pt); */

runtime_stat_t ret;

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
  /*   pt[m->method_id].mesh_timings = Init_array(smooth_numbers_n, sizeof(timing_t)); */
  /*   memset(pt[m->method_id].mesh_timings, 0, smooth_numbers_n * sizeof(timing_t)); */
  /* } */
  /* if((pt[m->method_id].cao_timings == NULL) || (pt[m->method_id].n != s->nparticles)) { */
  /*   puts("Resetting cao timings."); */
  /*   Free_array(pt[m->method_id].cao_timings); */
  /*   pt[m->method_id].cao_timings = Init_array(CAO_MAX+1, sizeof(runtime_stat_t)); */
  /*   memset(pt[m->method_id].cao_timings, 0, (CAO_MAX+1) * sizeof(runtime_stat_t)); */
  /*   pt[m->method_id].n = s->nparticles; */
  /* } */

  /* int *mesh_of = bsearch(&(p->mesh), smooth_numbers, smooth_numbers_n, sizeof(int), compare_ints); */
  /* int mesh_id = mesh_of - smooth_numbers; */

  /* if(pt[m->method_id].mesh_timings[mesh_id].avg <= 0) { */
  /*   pt[m->method_id].mesh_timings[mesh_id] = time_mesh(m,s,p); */
  /*   /\* TUNE_TRACE(printf("mesh miss %d (id %d), time %e\n", p->mesh, mesh_id, pt[m->method_id].mesh_timings[mesh_id]);); *\/ */
  /* } */
  /* if(pt[m->method_id].cao_timings[p->cao].t_c.avg <= 0) { */
  /*   pt[m->method_id].cao_timings[p->cao] = time_cao(m,s,p); */
  /*   /\* TUNE_TRACE(printf("cao miss %d, time %e\n", p->cao, s->nparticles * pt[m->method_id].cao_timings[p->cao]);); *\/ */
  /* } */

 ret = time_full(m, s, p);

 /* ret = pt[m->method_id].cao_timings[p->cao]; */
 /* ret.t_g = pt[m->method_id].mesh_timings[mesh_id]; */

 /* ret.t.avg = ret.t_c.avg + ret.t_g.avg + ret.t_f.avg; */
 /* ret.t.min = ret.t_c.min + ret.t_g.min + ret.t_f.min; */
 /* ret.t.max = ret.t_c.max + ret.t_g.max + ret.t_f.max; */

 return  ret;
}

runtime_stat_t Tune( const method_t *m, system_t *s, parameters_t *p, FLOAT_TYPE precision ) {
  //Parameter iteraters, best parameter set
  parameters_t it, p_best;
  // Array to store forces
// @TODO: Make copy of reference forces.

  forces_t *f = Init_forces(s->nparticles);

  it.prefactor = 1.0;

  runtime_stat_t best_time;
  best_time.t.avg = DBL_MAX;
  best_time.t.min = DBL_MAX;
  runtime_stat_t time;

  FLOAT_TYPE error = -1.0;
  FLOAT_TYPE V = pow( s->length, 3);
  int cao_start;
  int mesh_it_min=0, mesh_it_max=smooth_numbers_n, mesh_it;
  int success_once = 0;
  int cao_min, cao_max, cao_limit, cao_last;

  runtime_stat_t ret;
  ret.t.avg = -1;

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
    runtime_stat_t min_time;
    min_time = get_timing(m, s, &it);

    if( min_time.t.min > best_time.t.min) {
      TUNE_TRACE(printf("Best possible time for (%d %d) = %e slower than best %e\n", it.mesh, it.cao, min_time.t.avg, best_time.t.avg););
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
     
      TUNE_TRACE(printf("\n Timing mesh %d cao %d rcut %e time %e prec %e alpha %e\n", it.mesh, it.cao, it.rcut, time.t.avg, error, it.alpha ););
	
      if( time.t.min < best_time.t.min ) {
	p_best = it;
	p_best.precision = error;
	best_time = time;
	ret = best_time;
	TUNE_TRACE(printf("New best. mesh %d cao %d rcut %e time %e prec %e alpha %e\n", p_best.mesh, p_best.cao, p_best.rcut, time.t.avg, error, p_best.alpha );)         ;
      }
    }
  }
  Free_forces(f);

  if( success_once == 0 )
    return ret;

  *p = p_best;
  p->ip = p->cao - 1;
  p->cao3 = p->cao * p->cao * p->cao;

  TUNE_TRACE(printf("Using mesh %d cao %d rcut %e alpha %e with precision %e time %e\n", p->mesh, p->cao, p->rcut, p->alpha, p->precision, best_time.t.avg);)
  return ret;
}


