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
