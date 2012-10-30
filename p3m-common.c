#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>

#include "p3m-common.h"

#include "common.h"
#include "interpol.h"

int P3M_BRILLOUIN = 1;
int P3M_BRILLOUIN_TUNING = 1;

#define FREE_TRACE(A) 


FLOAT_TYPE sinc(FLOAT_TYPE d)
{
  double PId = PI*d;
  return (d == 0.0) ? 1.0 : SIN(PId)/PId;
}

FLOAT_TYPE analytic_cotangent_sum(int n, FLOAT_TYPE mesh_i, int cao)
{
    FLOAT_TYPE c, res=0.0;
    c = SQR(COS(PI*mesh_i*(FLOAT_TYPE)n));

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
        dn   -= ROUND(dn/dMesh)*dMesh;
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
        d->nshift[i] = i - ROUND(i/dMesh)*dMesh;

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
      printf("max %d\n", max);
      d->cf[1] = NULL;
      d->ca_ind[1] = NULL;
      
      for (i = 0; i < max; i++) {
	d->cf[i] = Init_array( p->cao3 * s->nparticles, sizeof(FLOAT_TYPE));
	d->ca_ind[i] = Init_array( 3*s->nparticles, sizeof(int));
      }
	
      d->inter = Init_interpolation( p->ip, m->flags & METHOD_FLAG_ad );
	
    }
    else {
      d->cf[0] = NULL;
      d->ca_ind[0] = NULL;
      d->cf[1] = NULL;
      d->ca_ind[1] = NULL;
      d->inter = NULL;
    }

    if ( m->flags & METHOD_FLAG_G_hat) {
        d->G_hat = Init_array(mesh3, sizeof(FLOAT_TYPE));
        m->Influence_function( s, p, d );   
    }
    else
        d->G_hat = NULL;    


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
        FFTW_FREE(d->G_hat);

    FREE_TRACE(puts("Free qmesh.");)
    if (d->Qmesh != NULL)
        FFTW_FREE(d->Qmesh);

    FREE_TRACE(puts("Free Fmesh.");)
    Free_vector_array(d->Fmesh);

    FREE_TRACE(puts("Free dshift.");)
    if (d->nshift != NULL)
        FFTW_FREE(d->nshift);

    if (d->Dn != NULL)
        FFTW_FREE(d->Dn);

    for (i=0;i<2;i++) {
        if (d->dQdx[i] != NULL)
            FFTW_FREE(d->dQdx[i]);
        if (d->dQdy[i] != NULL)
            FFTW_FREE(d->dQdy[i]);
        if (d->dQdz[i] != NULL)
            FFTW_FREE(d->dQdz[i]);
    }

    if(d->inter != NULL ) {
      Free_interpolation(d->inter);
      d->inter = NULL;
    }

    for (i=0;i<2;i++) {
        if (d->cf[i] != NULL)
            FFTW_FREE(d->cf[i]);
        if (d->ca_ind[i] != NULL)
            FFTW_FREE(d->ca_ind[i]);
    }

    for(i=0; i<d->forward_plans; i++) {
      FFTW_DESTROY_PLAN(d->forward_plan[i]);
    }

    for(i=0; i<d->backward_plans; i++) {
      FFTW_DESTROY_PLAN(d->backward_plan[i]);
    }

    FFTW_FREE(d);

}

FLOAT_TYPE C_ewald(int nx, int ny, int nz, system_t *s, parameters_t *p) {
  int mx, my, mz;
  int nmx, nmy, nmz;
  FLOAT_TYPE km2;
  FLOAT_TYPE ret = 0.0;

  for (mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++) {
    nmx = nx + p->mesh*mx;
    for (my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++) {
      nmy = ny + p->mesh*my;
      for (mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++) {
	nmz = nz + p->mesh*mz;

	km2 = SQR(2.0*PI/s->length) * ( SQR ( nmx ) + SQR ( nmy ) + SQR ( nmz ) );
	ret += EXP(- 2.0 * km2 / ( 4.0 * SQR(p->alpha)) ) / km2;
      }
    }
  }
  return 16.0 * SQR(PI) * ret;
}

#define NTRANS(N) (N<0) ? (N + d->mesh) : N

FLOAT_TYPE Generic_error_estimate(R3_to_R A, R3_to_R B, R3_to_R C, system_t *s, parameters_t *p, data_t *d) {
  // The Hockney-Eastwood pair-error functional.
  FLOAT_TYPE Q_HE = 0.0, Q_opt=0.0;
  // Linear index for G, this breaks notation, but G is calculated anyway, so this is convinient.
  int ind = 0;
  // Convinience variable to hold the current value of the influence function.
  FLOAT_TYPE G_hat = 0.0;
  FLOAT_TYPE V = s->length * SQR(s->length);
  FLOAT_TYPE a,b,c;

  int nx, ny, nz;

  for (nx=-d->mesh/2; nx<d->mesh/2; nx++) {
    for (ny=-d->mesh/2; ny<d->mesh/2; ny++) {
      for (nz=-d->mesh/2; nz<d->mesh/2; nz++) {
	if((nx!=0) && (ny!=0) && (nz!=0)) {
	  ind = r_ind(NTRANS(nx), NTRANS(ny), NTRANS(nz));
	  G_hat = d->G_hat[ind];

	  a = A(nx,ny,nz,s,p);
	  b = B(nx,ny,nz,s,p);
	  c = C(nx,ny,nz,s,p);

	  Q_HE += a * SQR(G_hat) - 2.0 * b * G_hat + c;
	  Q_opt += c - SQR( b ) / a;
	  printf("A\t%e\tB\t%e\tC\t%e\n", FLOAT_CAST a , FLOAT_CAST b, FLOAT_CAST c);
	  printf("B/A\t%e\tG_hat\t%e\n", FLOAT_CAST (b / a), FLOAT_CAST G_hat);
	  printf("Q_HE\t%lf\tQ_opt\t%lf\n", FLOAT_CAST Q_HE, FLOAT_CAST Q_opt);
	}
      }
    }
  }
  printf("Final Q_HE\t%lf\tQ_opt\t%lf\n", FLOAT_CAST Q_HE, FLOAT_CAST Q_opt);
  printf("dF_opt = %e\n", FLOAT_CAST (s->q2* SQRT( FLOAT_ABS(Q_opt) / (FLOAT_TYPE)s->nparticles) / V));

  return  s->q2* SQRT( FLOAT_ABS(Q_HE) / (FLOAT_TYPE)s->nparticles) / V;

}
