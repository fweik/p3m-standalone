#include <math.h>
#include <stdlib.h>

#include "p3m-common.h"

FLOAT_TYPE sinc(FLOAT_TYPE d)
{
#define epsi 0.1

#define c2 -0.1666666666667e-0
#define c4  0.8333333333333e-2
#define c6 -0.1984126984127e-3
#define c8  0.2755731922399e-5

  double PId = PI*d, PId2;

  if (fabs(d)>epsi)
    return sin(PId)/PId;
  else {
    PId2 = SQR(PId);
    return 1.0 + PId2*(c2+PId2*(c4+PId2*(c6+PId2*c8)));
  }
  return 1.0;
}

FLOAT_TYPE analytic_cotangent_sum(int n, FLOAT_TYPE mesh_i, int cao)
{
  FLOAT_TYPE c, res=0.0;
  c = SQR(cos(PI*mesh_i*(FLOAT_TYPE)n));

  switch (cao) {
  case 1 : { 
    res = 1; 
    break; }
  case 2 : { 
    res = (1.0+c*2.0)/3.0; 
    break; }
  case 3 : { 
    res = (2.0+c*(11.0+c*2.0))/15.0; 
    break; }
  case 4 : { 
    res = (17.0+c*(180.0+c*(114.0+c*4.0)))/315.0; 
    break; }
  case 5 : { 
    res = (62.0+c*(1072.0+c*(1452.0+c*(247.0+c*2.0))))/2835.0; 
    break; }
  case 6 : { 
    res = (1382.0+c*(35396.0+c*(83021.0+c*(34096.0+c*(2026.0+c*4.0)))))/155925.0; 
    break; }
  case 7 : { 
    res = (21844.0+c*(776661.0+c*(2801040.0+c*(2123860.0+c*(349500.0+c*(8166.0+c*4.0))))))/6081075.0; 
    break; }
  }
  
  return res;
}


void Differenzenoperator_berechnen(p3m_parameters_t *p, p3m_data_t *d)
{
  /* 
     Die Routine berechnet den fourieretransformierten 
     Differentialoperator auf der Ebene der n, nicht der k,
     d.h. der Faktor  i*2*PI/L fehlt hier!
  */
  
  int    i;
  FLOAT_TYPE dMesh=(FLOAT_TYPE)d->mesh;
  FLOAT_TYPE dn;

  d->Dn = Init_array(d->mesh, sizeof(FLOAT_TYPE));  

  for (i=0; i<d->mesh; i++)
    {
      dn    = (FLOAT_TYPE)i; 
      dn   -= round(dn/dMesh)*dMesh;
      d->Dn[i] = dn;
    }

  d->Dn[d->mesh/2] = 0.0;  

}

void nshift_ausrechnen(p3m_data_t *d)
{
  /* Verschiebt die Meshpunkte um Mesh/2 */

  int    i;
  FLOAT_TYPE dMesh=(FLOAT_TYPE)d->mesh;

  d->nshift = Init_array(d->mesh, sizeof(FLOAT_TYPE));  

  for (i=0; i<d->mesh; i++) 
    d->nshift[i] = i - round(i/dMesh)*dMesh; 
  
}

void Init_data(method_t *m, p3m_parameters_t *p, p3m_data_t *d) {
  
};

void Free_data(p3m_data_t *d) {
  int i;
  if(d->G_hat != NULL)
    free(d->G_hat);
  if(d->Qmesh != NULL)
    free(d->Qmesh);

  Free_vector_array(&d->Fmesh);
  
  if(d->nshift != NULL)
    free(d->nshift);
  if(d->Dn != NULL)
    free(d->Dn);
  for(i=0;i<2;i++) {
    if(d->dQdx[i] != NULL)
      free(d->dQdx[i]);  
    if(d->dQdy[i] != NULL)
      free(d->dQdy[i]);  
    if(d->dQdz[i] != NULL)
      free(d->dQdz[i]);  
  }
}