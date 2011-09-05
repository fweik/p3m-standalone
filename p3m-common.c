#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>

#include "p3m-common.h"

#include "common.h"
#include "interpol.h"


#define FREE_TRACE(A) 

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
    case 1 :
        res = 1;
        break;
    case 2 :
        res = (1.0+c*2.0)/3.0;
        break;
    case 3 :
        res = (2.0+c*(11.0+c*2.0))/15.0;
        break;
    case 4 :
        res = (17.0+c*(180.0+c*(114.0+c*4.0)))/315.0;
        break;
    case 5 :
        res = (62.0+c*(1072.0+c*(1452.0+c*(247.0+c*2.0))))/2835.0;
        break;
    case 6 :
        res = (1382.0+c*(35396.0+c*(83021.0+c*(34096.0+c*(2026.0+c*4.0)))))/155925.0;
        break;
    case 7 :
        res = (21844.0+c*(776661.0+c*(2801040.0+c*(2123860.0+c*(349500.0+c*(8166.0+c*4.0))))))/6081075.0;
        break;
    }

    return res;
}


void Init_differential_operator(data_t *d)
{
    /*
       Die Routine berechnet den fourieretransformierten
       Differentialoperator auf er Ebene der n, nicht der k,
       d.h. der Faktor  i*2*PI/L fehlt hier!
    */

    int    i;
    FLOAT_TYPE dMesh=(FLOAT_TYPE)d->mesh;
    FLOAT_TYPE dn;

    for (i=0; i<d->mesh; i++)
    {
        dn    = (FLOAT_TYPE)i;
        dn   -= round(dn/dMesh)*dMesh;
        d->Dn[i] = dn;
    }

    d->Dn[d->mesh/2] = 0.0;

}

void Init_nshift(data_t *d)
{
    /* Verschiebt die Meshpunkte um Mesh/2 */

    int    i;
    FLOAT_TYPE dMesh=(FLOAT_TYPE)d->mesh;

    for (i=0; i<d->mesh; i++)
        d->nshift[i] = i - round(i/dMesh)*dMesh;

}

data_t *Init_data(const method_t *m, system_t *s, parameters_t *p) {
    int mesh3 = p->mesh*p->mesh*p->mesh;
    data_t *d = Init_array(1, sizeof(data_t));

    d->mesh = p->mesh;
    
    if ( m->flags & METHOD_FLAG_Qmesh)
        d->Qmesh = Init_array(2*mesh3, sizeof(FLOAT_TYPE));
    else
        d->Qmesh = NULL;

    if ( m->flags & METHOD_FLAG_ik ) {
        d->Fmesh = Init_vector_array(2*mesh3);
        d->Dn = Init_array(d->mesh, sizeof(FLOAT_TYPE));
        Init_differential_operator(d);
    }
    else {
        d->Fmesh = NULL;
        d->Dn = NULL;
    }

    d->nshift = NULL;

    if ( m->flags & METHOD_FLAG_nshift ) {
      d->nshift = Init_array(d->mesh, sizeof(FLOAT_TYPE));
      Init_nshift(d);
    }

    if ( m->flags & METHOD_FLAG_G_hat) {
        d->G_hat = Init_array(mesh3, sizeof(FLOAT_TYPE));
        m->Influence_function( s, p, d );   
    }
    else
        d->G_hat = NULL;    

    d->dQdx[0] = d->dQdy[0] = d->dQdz[0] = NULL;
    d->dQdx[1] = d->dQdy[1] = d->dQdz[1] = NULL;

    if ( m->flags & METHOD_FLAG_ad ) {
      int i;
      int max = ( m->flags & METHOD_FLAG_interlaced) ? 2 : 1;

        for (i = 0; i < max; i++) {
            d->dQdx[i] = Init_array( s->nparticles*p->cao3, sizeof(FLOAT_TYPE) );
            d->dQdy[i] = Init_array( s->nparticles*p->cao3, sizeof(FLOAT_TYPE) );
            d->dQdz[i] = Init_array( s->nparticles*p->cao3, sizeof(FLOAT_TYPE) );
        }
    }

    if ( m->flags & METHOD_FLAG_ca ) {
      int i;
      int max = ( m->flags & METHOD_FLAG_interlaced ) ? 2 : 1;
      
      d->cf[1] = NULL;
      d->ca_ind[1] = NULL;
      
      for (i = 0; i < max; i++) {
	d->cf[i] = Init_array( p->cao3 * s->nparticles, sizeof(FLOAT_TYPE));
	d->ca_ind[i] = Init_array( 3*s->nparticles, sizeof(int));
      }
	
      Init_interpolation( p->ip, d );
	
    }
    else {
      d->cf[0] = NULL;
      d->ca_ind[0] = NULL;
      d->cf[1] = NULL;
      d->ca_ind[1] = NULL;
      d->LadInt = NULL;
      d->LadInt_ = NULL;
    }

    d->forward_plans = 0;
    d->backward_plans = 0;

    return d;
}

void Free_data(data_t *d) {
    int i;

    if( d == NULL )
      return;

    FREE_TRACE(puts("Free_data(); Free ghat.");)
    if (d->G_hat != NULL)
        fftw_free(d->G_hat);

    FREE_TRACE(puts("Free qmesh.");)
    if (d->Qmesh != NULL)
        fftw_free(d->Qmesh);

    FREE_TRACE(puts("Free Fmesh.");)
    Free_vector_array(d->Fmesh);

    FREE_TRACE(puts("Free dshift.");)
    if (d->nshift != NULL)
        fftw_free(d->nshift);

    if (d->Dn != NULL)
        fftw_free(d->Dn);

    for (i=0;i<2;i++) {
        if (d->dQdx[i] != NULL)
            fftw_free(d->dQdx[i]);
        if (d->dQdy[i] != NULL)
            fftw_free(d->dQdy[i]);
        if (d->dQdz[i] != NULL)
            fftw_free(d->dQdz[i]);
    }

    if(d->LadInt != NULL ) {
      for(i=0;d->LadInt[i] != NULL; i++)
	fftw_free(d->LadInt[i]);
      fftw_free(d->LadInt);
    }

    if(d->LadInt_ != NULL ) {
      for(i=0;d->LadInt_[i] != NULL; i++)
	fftw_free(d->LadInt_[i]);
      fftw_free(d->LadInt_);
    }

    for (i=0;i<2;i++) {
        if (d->cf[i] != NULL)
            fftw_free(d->cf[i]);
        if (d->ca_ind[i] != NULL)
            fftw_free(d->ca_ind[i]);
    }

    for(i=0; i<d->forward_plans; i++) {
      fftw_destroy_plan(d->forward_plan[i]);
    }

    for(i=0; i<d->backward_plans; i++) {
      fftw_destroy_plan(d->backward_plan[i]);
    }

    fftw_free(d);

}
