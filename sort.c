#include "sort.h"

void sort_particles_r(system_t *s, int i, int j);

// Swap particles i and j in s
inline void swap_particles( system_t *s, int i, int j) {
  FLOAT_TYPE b;
  
  b = s->p->x[j];
  s->p->x[j] = s->p->x[i];
  s->p->x[i] = b;

  b = s->p->y[j];
  s->p->y[j] = s->p->y[i];
  s->p->y[i] = b;

  b = s->p->z[j];
  s->p->z[j] = s->p->z[i];
  s->p->z[i] = b;

  b = s->q[j];
  s->q[j] = s->q[i];
  s->q[i] = b;
}

// Sort system particles by x coordinate, using quicksort ( O(n*log(n)) average )
void sort_particles(system_t *s) {
  sort_particles_r( s, 0, s->nparticles );
}

void sort_particles_r(system_t *s, int i, int j) {
  int pivot = (i+j)/2;
  int m,n;
  FLOAT_TYPE key;
  FLOAT_TYPE *x = s->p->x;

  if( i < j ) {
    swap_particles( s, i, pivot );
    key = x[pivot];

    m=i+1;
    n = j;
    while(m<=n) {
      while((m<=j) && (x[m] <= key))
	m++;
      while((n>=i) && (x[n] > key))
	j--;
      if(i<j)
	swap_particles(s, m, n);

      swap_particles(s, i, n);
      sort_particles_r( s, i, n - 1);
      sort_particles_r( s, n+1, j );
    }
  }
}
