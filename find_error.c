#include <stdio.h>
#include <math.h>

#include "q_ik_double.h"

double p3m_find_alpha(double q, int mesh_id, int cao_id) {
  const double *p = Q_ik[mesh_id][cao_id];
  double h = alphaLmax / (points-1);

  if(p[0] > q)
    return -1;

  int a = 0, b = points-1;
  int c;
  
  while((b-a) > 1) {
    c = (a+b)/2;

    if(fabs(p[c] - q) < 1e-15) {
      return c*h;
    }

    if(p[c] < q)
      a = c;
    else
      b = c;
  }
  return a*h;
}

double p3m_find_error(double alphaL, int mesh_id, int cao_id) {
  double h = alphaLmax / (points-1);
  double d;

  d = alphaL/h;
  int l = (int)floor(d);
  d -= l;

  return (1.-d)*Q_ik[mesh_id][cao_id][l] + d*Q_ik[mesh_id][cao_id][l];
}  
