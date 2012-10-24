#include "sort.h"

void sort_particles_r(system_t *s, int m, int n);

// Swap particles i and j in s
inline static void swap_particles( system_t *s, int i, int j) {
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
  sort_particles_r( s, 0, s->nparticles - 1 );
}

void sort_particles_r(system_t *s, int m, int n) {
  int k = (m+n)/2;
  int i,j;
  FLOAT_TYPE key;
  FLOAT_TYPE *x = s->p->x;

  if( m < n ) {
    swap_particles( s, m, k );
    key = x[m];

    i = m+1;
    j = n;
    while( i <= j) {
      while(( i <= n ) && (x[i] <= key))
	i++;
      while(( j >= m ) && (x[j] > key))
	j--;

      if( i < j )
	swap_particles(s, i, j);
    }
      
    swap_particles( s, m, j );
      
    sort_particles_r( s, m, j - 1);
    sort_particles_r( s, j+1 , n );
    
  }
}
