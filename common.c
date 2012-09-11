#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <fftw3.h>

#ifdef DETAILED_TIMINGS
#include <mpi.h>
#endif

#include "types.h"
#include "common.h"
#include "realpart.h"
#include "ewald.h"
#include "p3m-common.h"

#ifdef __detailed_timings
double t_charge_assignment[4];
double t_force_assignment[4];
double t_real_part[4];
double t_k_part[4];
double t_convolution[4];
double t_fft[4];
double timer;
#endif


void *Init_array(int size, size_t field_size) {
    void *a;

    assert(size >= 0);
    assert(field_size > 0);

    a = FFTW_MALLOC(size * field_size);

    assert(a != NULL);
    return a;
}

vector_array_t *Init_vector_array(int n) {
    int i;
    vector_array_t *v;

    v = Init_array( 1, sizeof(vector_array_t));

    assert(n >= 0 );

    v->fields = Init_array( 3, sizeof(FLOAT_TYPE *));


    assert( v->fields != NULL );

    for (i=0;i<3;i++) {
        v->fields[i] = Init_array( n, sizeof(FLOAT_TYPE));
    }

    v->x = v->fields[0];
    v->y = v->fields[1];
    v->z = v->fields[2];

    return v;
}

system_t *Init_system(int n) {
    system_t *s;

    assert( n > 0 );

    s = Init_array( 1, sizeof(system_t));

    s-> nparticles = n;

    s-> p = Init_vector_array(s->nparticles);

    s->reference = Init_forces(s->nparticles);

    s->q = Init_array(s->nparticles, sizeof(FLOAT_TYPE));

    return s;
}

forces_t *Init_forces(int n) {
    forces_t *f;

    f = Init_array( 1, sizeof(forces_t));

    f->f = Init_vector_array(n);
    f->f_k = Init_vector_array(n);
    f->f_r = Init_vector_array(n);
    
    return f;
}

void Free_system(system_t *s) {
  if( s == NULL)
    return;
  
  Free_vector_array(s->p);
  
  Free_forces(s->reference);
  
  FFTW_FREE(s->q);
  
  FFTW_FREE(s);
  
}

void Free_forces(forces_t *f) {
  if( f == NULL )
    return;
  
  Free_vector_array(f->f);
  Free_vector_array(f->f_k);
  Free_vector_array(f->f_r);
 
  FFTW_FREE(f);
}

void Free_vector_array(vector_array_t *v) {
    int i;

    if ( v != NULL ) {
        for (i=0;i<3;i++) {
            if (v->fields[i] != NULL) {
                FFTW_FREE(v->fields[i]);
                v->fields[i] = NULL;
            }
        }

        FFTW_FREE(v->fields);

        v->x = v->y = v->z = NULL;
        v->fields = NULL;

	FFTW_FREE(v);
    }
}

void Calculate_forces ( const method_t *m, system_t *s, parameters_t *p, data_t *d, forces_t *f ) {

    int i, j;

#ifdef DETAILED_TIMINGS
    FLOAT_TYPE t;
#endif

    for ( i=0; i<3; i++ ) {
        memset ( f->f->fields[i]  , 0, s->nparticles*sizeof ( FLOAT_TYPE ) );
        memset ( f->f_k->fields[i], 0, s->nparticles*sizeof ( FLOAT_TYPE ) );
        memset ( f->f_r->fields[i], 0, s->nparticles*sizeof ( FLOAT_TYPE ) );
    }

#ifdef DETAILED_TIMINGS
    if(__detailed_timings)
      t = MPI_Wtime();
#endif
    //Realpart_neighborlist ( s, p, d, f );
#ifdef DETAILED_TIMINGS
    if(__detailed_timings)
      printf(" %e", MPI_Wtime() - t);
#endif


    Realteil( s, p, f );
    //Realpart_neighborlist( s, p, d, f );

    //  Dipol(s, p);

    m->Kspace_force ( s, p, d, f );

    //#pragma omp parallel for private( i )
    for ( j=0; j < 3; j++ ) {
        for ( i=0; i<s->nparticles; i++ ) {
            f->f->fields[j][i] += f->f_k->fields[j][i] + f->f_r->fields[j][i];
        }
    }
}

FLOAT_TYPE Calculate_reference_forces ( system_t *s, parameters_t *p ) {

    data_t *d;

    parameters_t op = *p;

    forces_t *f = Init_forces ( s->nparticles );

    op.rcut = 0.49 * s->length;

    op.alpha = Ewald_compute_optimal_alpha ( s, &op );

    printf("Reference Forces: alpha %lf, r_cut %lf\n", op.alpha, op.rcut);

    d = method_ewald.Init ( s, &op );

    //Init_neighborlist( s, &op, d );
     
    method_ewald.Influence_function( s, &op, d );

    s->energy = 0.0;
    
    Calculate_forces ( &method_ewald, s, &op, d, s->reference );

    //Free_neighborlist(d);
    Free_data(d);
    Free_forces(f);

    return method_ewald.Error( s, &op );
}

FLOAT_TYPE distance( system_t *s, int i, int j ) {
  int k;
  FLOAT_TYPE ret = 0.0;
  for(k=0;k<3;k++) {
    ret += SQR(s->p->fields[k][i] - s->p->fields[k][j]);
  }
  return SQRT(ret);
}

FLOAT_TYPE Min_distance( system_t *s ) {
  int i,j;
  FLOAT_TYPE min =2.0*s->length, d;
  

  for(i=0;i<s->nparticles;i++) {
    for(j=0;j<i;j++) {
      d = distance( s, i, j);
      min = (d < min) ? d : min;
    }
  }
  return min;
}
