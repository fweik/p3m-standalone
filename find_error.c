#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "q.h"

static int
compare_ints (const void *a, const void *b)
{
  return (*(int *)(a) - *(int *)(b));
}

double p3m_find_alpha(double q, int mesh_id, int cao_id) {
  const double *p = Q_ik[mesh_id][cao_id];
  double h = alphaLmax / (n_points-1);

  if(p[0] > q)
    return -1;

  int a = 0, b = n_points-1;
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

double p3m_find_q_ik(double alphaL, int mesh_id, int cao_id) {
  double h = alphaLmax / (n_points-1);
  double d;

  d = alphaL/h;
  int l = (int)floor(d);
  d -= l;

  return (1.-d)*Q_ik[mesh_id][cao_id][l] + d*Q_ik[mesh_id][cao_id][l+1];
}  

double p3m_find_q_ik_i(double alphaL, int mesh_id, int cao_id) {
  double h = alphaLmax / (n_points-1);
  double d;

  d = alphaL/h;
  int l = (int)floor(d);
  d -= l;

  return (1.-d)*Q_ik_i[mesh_id][cao_id][l] + d*Q_ik_i[mesh_id][cao_id][l+1];
}  

double p3m_find_q_ad(double alphaL, int mesh_id, int cao_id) {
  double h = alphaLmax / (n_points-1);
  double d;

  d = alphaL/h;
  int l = (int)floor(d);
  d -= l;

  return (1.-d)*Q_ad[mesh_id][cao_id][l] + d*Q_ad[mesh_id][cao_id][l+1];
}  

double p3m_find_q_ad_i(double alphaL, int mesh_id, int cao_id) {
  double h = alphaLmax / (n_points-1);
  double d;

  d = alphaL/h;
  int l = (int)floor(d);
  d -= l;

  return (1.-d)*Q_ad_i[mesh_id][cao_id][l] + d*Q_ad_i[mesh_id][cao_id][l+1];
}  



double p3m_find_error(double alphaL, int mesh, int cao, int method) {
    int *mesh_id, *cao_id;
    // Check whether value pair is tabulated. 
  
    mesh_id = bsearch(&mesh, meshes, n_meshes, sizeof(int), compare_ints);
    cao_id =  bsearch(&cao, caos, n_caos, sizeof(int), compare_ints);

    // Parameter set not tabulated
    if((mesh_id == NULL) || (cao_id == NULL) || (alphaL > alphaLmax))
      return -1.;
    else
      switch(method) {
      case 0:
	return p3m_find_q_ik(alphaL, mesh_id - meshes, cao_id - caos);
      case 1:
	return p3m_find_q_ik_i(alphaL, mesh_id - meshes, cao_id - caos);
      case 2: 
	return p3m_find_q_ad(alphaL, mesh_id - meshes, cao_id - caos); 
      case 3:
	return p3m_find_q_ad_i(alphaL, mesh_id - meshes, cao_id - caos); 
      default:
	return -1;
      }
}
