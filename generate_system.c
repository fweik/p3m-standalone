#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

#include "generate_system.h"

#include "common.h"

#include "math.h"

#include "stdlib.h"

system_t *generate_seperated_dipole(int size, FLOAT_TYPE box) {
  int id=0,i,j;
  system_t *s;
  FLOAT_TYPE d = 1.0;
  int n_dipoles = size / 2;
  FLOAT_TYPE u, theta;
  FLOAT_TYPE x[3], y[3];
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

  s = Init_system(2*n_dipoles);
  s->length = box;
  s->q2 = 0.0;

  /* Generate randome position for the first particle of the dipole x,
     then choose random direction and place second particle at x + d*y. Second
     particle could be outside the box and has to be folded. */

  for(i=0; i < n_dipoles; i++) {
    for(j=0;j<3;j++) {
      x[j] = box * gsl_rng_uniform(rng);
    }
    u = -1 + 2*gsl_rng_uniform(rng);
    theta = 2*PI*gsl_rng_uniform(rng);

    y[0] = SQRT(1 - SQR(u)) * COS(theta);
    y[1] = SQRT(1 - SQR(u)) * SIN(theta);
    y[2] = u;

    for(j=0;j<3;j++) {
      s->p->fields[j][id  ] = x[j];
      s->p->fields[j][id+1] = x[j] + d*y[j];
      while(s->p->fields[j][id+1] >= box)
	s->p->fields[j][id+1] -= box;
      while(s->p->fields[j][id+1] < 0.0)
	s->p->fields[j][id+1] += box;
    }
    s->q[id  ] = +1.0;
    s->q[id+1] = -1.0;
    s->q2 += SQR(s->q[id]) + SQR(s->q[id+1]);

    id+=2;  
  }

  return s;
}

system_t *generate_madelung(int size, FLOAT_TYPE box) {
  int id=0,i,j,k;
  system_t *s;
  int per_row = ceil( pow(size, 0.33333) );
  FLOAT_TYPE a = box / per_row;
  FLOAT_TYPE off = 0.5*a;

  size = per_row * per_row * per_row; 

  printf("generate_madelung: box %lf, size %d, per_row %d, a %lf, off %lf\n", box, size, per_row, a, off);

  s = Init_system(size);
  s->length = box;
  s->q2 = 0.0;

  for(i=0; i<per_row; i++)
    for(j=0; j<per_row; j++)
      for(k=0; k<per_row; k++) {
	s->p->fields[0][id] = a*i + off;
	s->p->fields[1][id] = a*j + off;
	s->p->fields[2][id] = a*k + off;

	s->q[id] = 1.0 - 2.0 * ((i+j+k) % 2);
	s->q2 += s->q[id] * s->q[id];
	id++;
      }
  return s;
}

system_t *generate_inner_box(int size, FLOAT_TYPE box) {
  FLOAT_TYPE lower_left = 0.25 * box;
  puts("Generating inner box system.");
  int i,j;
  system_t *s;
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

  s = Init_system(size);
  s->length = box;
  s->q2 = 0.0;

  for(i=0;i<size;i++) {
    for(j=0;j<3;j++) {
      s->p->fields[j][i] = lower_left + 0.5 * box*gsl_rng_uniform(rng);
    }
    s->q[i] = 1.0 - 2.0 * (i%2);
    s->q2 += s->q[i] * s->q[i];
  }
  return s;
}

system_t *generate_random_system(int size, FLOAT_TYPE box, FLOAT_TYPE max_charge ) {
  int i,j;
  system_t *s;
  FLOAT_TYPE q2 = 0.0;
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  
  s = Init_system(size);
  s->length = box;
  s->q2 = 0.0;

  for(i=0;i<size;i++) {
    for(j=0;j<3;j++) {
      s->p->fields[j][i] = box*gsl_rng_uniform(rng);
    }
    s->q[i] = (1.0 - 2.0 * (i%2)) * max_charge;
    q2 += s->q[i] * s->q[i];
  }
  s->q2 = q2;

  gsl_rng_free (rng);

  return s;
}

system_t *generate_system( int form_factor, int size, FLOAT_TYPE box, FLOAT_TYPE max_charge ) {
  srand48(42);
  switch(form_factor) {
  case FORM_FACTOR_RANDOM: 
    return generate_random_system( size, box, max_charge );
    break;
  case FORM_FACTOR_INNER_BOX:
    return generate_inner_box(size, box);
    break;
  case FORM_FACTOR_MADELUNG:
    return generate_madelung(size,box);
    break;
  case FORM_FACTOR_SEPERATED_DIPOLE:
    return generate_seperated_dipole(size,box);
    break;
  default:
    fprintf( stderr, "Warning, form factor not known.\n");
  }
  return (void *)NULL;
}
