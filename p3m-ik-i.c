
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

#include "types.h"
#include "common.h"
#include "p3m-common.h"
#include "charge-assign.h"
#include "p3m-ik-i.h"

#include "realpart.h"

#ifdef __detailed_timings
#include <mpi.h>
#endif


// declaration of the method

const method_t method_p3m_ik_i = { METHOD_P3M_ik_i, "P3M with ik differentiation, interlaced.",
                                   METHOD_FLAG_P3M | METHOD_FLAG_ik | METHOD_FLAG_interlaced,
                                   &Init_ik_i, &Influence_ik_i, &P3M_ik_i, &Error_ik_i, &p3m_k_space_error_ik_i,
                                 };

// Forward declaration of local functions

static void forward_fft( data_t * );
static void backward_fft( data_t * );

inline void forward_fft ( data_t *d ) {
    FFTW_EXECUTE ( d->forward_plan[0] );
}

inline void backward_fft ( data_t *d ) {
    int i;
    for ( i=0;i<3;i++ )
        FFTW_EXECUTE ( d->backward_plan[i] );
}

data_t *Init_ik_i ( system_t *s, parameters_t *p ) {
    int l;
    int mesh = p->mesh;

    data_t *d = Init_data ( &method_p3m_ik_i, s, p );

    d->forward_plans = 1;
    d->backward_plans = 3;

    d->forward_plan[0] = FFTW_PLAN_DFT_3D ( mesh, mesh, mesh, ( FFTW_COMPLEX * ) d->Qmesh, ( FFTW_COMPLEX * ) d->Qmesh, FFTW_FORWARD, FFTW_PATIENT );

    for ( l=0;l<3;l++ ) {
        d->backward_plan[l] = FFTW_PLAN_DFT_3D ( mesh, mesh, mesh, ( FFTW_COMPLEX * ) ( d->Fmesh->fields[l] ), ( FFTW_COMPLEX * ) ( d->Fmesh->fields[l] ), FFTW_BACKWARD, FFTW_PATIENT );
    }
    return d;
}


void Aliasing_sums_ik_i( system_t *s, parameters_t *p, data_t *d, int NX, int NY, int NZ,
                         FLOAT_TYPE *Zaehler, FLOAT_TYPE *Nenner1, FLOAT_TYPE *Nenner2)  {
  FLOAT_TYPE S1,S2,S3;
  FLOAT_TYPE fak1,fak2,zwi;
  int    MX,MY,MZ;
  FLOAT_TYPE NMX,NMY,NMZ;
  FLOAT_TYPE NM2;
  
  FLOAT_TYPE expo, TE;
  int Mesh = d->mesh;
  
  FLOAT_TYPE Leni = 1.0/s->length;
  fak1 = 1.0/(FLOAT_TYPE)Mesh;
  fak2 = SQR(PI/p->alpha);
  
  Zaehler[0] = Zaehler[1] = Zaehler[2] = *Nenner1 = *Nenner2 = 0.0;
  
  for ( MX = -P3M_BRILLOUIN; MX <= P3M_BRILLOUIN; MX++ ) {
    NMX = d->nshift[NX] + Mesh*MX;
    S1  = my_power( sinc(fak1*NMX), 2*p->cao );
    for ( MY = -P3M_BRILLOUIN; MY <= P3M_BRILLOUIN; MY++ ) {
      NMY = d->nshift[NY] + Mesh*MY;
      S2   = S1 * my_power( sinc(fak1*NMY), 2*p->cao );
      for ( MZ = -P3M_BRILLOUIN; MZ <= P3M_BRILLOUIN; MZ++ ) {
	NMZ = d->nshift[NZ] + Mesh*MZ;
	S3  = S2*my_power( sinc(fak1*NMZ), 2*p->cao );
	
	NM2 = SQR ( NMX*Leni ) + SQR ( NMY*Leni ) + SQR ( NMZ*Leni );
	*Nenner1 += S3;
	
	expo = fak2*NM2;
	TE = EXP ( -expo );
	zwi  = S3 * TE/NM2;
	Zaehler[0] += NMX*zwi*Leni;
	Zaehler[1] += NMY*zwi*Leni;
	Zaehler[2] += NMZ*zwi*Leni;
	
	if (((MX+MY+MZ)%2)==0)                                       //even term
	  *Nenner2 += S3;
	else                                                //odd term: minus sign!
	  *Nenner2 -= S3;
      }
    }
  }
}



void Influence_ik_i( system_t *s, parameters_t *p, data_t *d )
{
    int    NX,NY,NZ;
    FLOAT_TYPE Dnx,Dny,Dnz;
    FLOAT_TYPE dMesh,dMeshi;
    FLOAT_TYPE Zaehler[3]={0.0,0.0,0.0},Nenner1=0.0, Nenner2=0.0;
    FLOAT_TYPE zwi;

    int ind = 0;
    int Mesh = p->mesh;
    FLOAT_TYPE Leni = 1.0/s->length;

    dMesh = (FLOAT_TYPE)Mesh;
    dMeshi= 1.0/dMesh;

    for (NX=0; NX<Mesh; NX++)
    {
        for (NY=0; NY<Mesh; NY++)
        {
            for (NZ=0; NZ<Mesh; NZ++)
            {
                ind = r_ind( NX, NY, NZ );

                if ((NX==0) && (NY==0) && (NZ==0))
                    d->G_hat[ind]=0.0;
                else if ((NX%(Mesh/2) == 0) && (NY%(Mesh/2) == 0) && (NZ%(Mesh/2) == 0))
                    d->G_hat[ind]=0.0;
                else
                {
                    Aliasing_sums_ik_i( s, p, d, NX, NY, NZ, Zaehler, &Nenner1,  &Nenner2);

                    Dnx = d->Dn[NX];
                    Dny = d->Dn[NY];
                    Dnz = d->Dn[NZ];

                    zwi  = Dnx*Zaehler[0]*Leni + Dny*Zaehler[1]*Leni + Dnz*Zaehler[2]*Leni;
                    zwi /= ( SQR(Dnx*Leni) + SQR(Dny*Leni) + SQR(Dnz*Leni) );
                    zwi /= 0.5*(SQR(Nenner1) + SQR(Nenner2));

                    d->G_hat[ind] = 2.0 * zwi / PI;
                }
            }
        }
    }
}


void P3M_ik_i( system_t *s, parameters_t *p, data_t *d, forces_t *f )
{
    /* Zaehlvariablen: */
    int i, j, k, l;
    /* Schnelles Modulo: */

    FLOAT_TYPE T1;

    FLOAT_TYPE Mesh = p->mesh;
    FLOAT_TYPE Leni = 1.0/s->length;

    FLOAT_TYPE dop;

    int c_index;

    /* Initialisieren von Qmesh */
    memset ( d->Qmesh, 0, 2*Mesh*Mesh*Mesh*sizeof ( FLOAT_TYPE ) );

  #ifdef __detailed_timings
  timer = MPI_Wtime();
  #endif


    /* chargeassignment */
    assign_charge( s, p, d, 0 );
    assign_charge( s, p, d, 1 );

  #ifdef __detailed_timings
  timer = MPI_Wtime() - timer;
  t_charge_assignment[1] = timer;
  timer = MPI_Wtime();
  #endif


    /* Durchfuehren der Fourier-Hin-Transformationen: */
    forward_fft(d);

  #ifdef __detailed_timings
    timer = MPI_Wtime() - timer;
    t_fft[1] = timer;
   timer = MPI_Wtime();
  #endif


    for (i=0; i<Mesh; i++)
        for (j=0; j<Mesh; j++)
            for (k=0; k<Mesh; k++)
            {
                c_index = c_ind(i,j,k);

                T1 = d->G_hat[r_ind(i,j,k)];
                d->Qmesh[c_index] *= T1;
                d->Qmesh[c_index+1] *= T1;

                for (l=0;l<3;l++) {
                    switch ( l ) {
                    case 0:
                        dop = d->Dn[i];
                        break;
                    case 1:
                        dop = d->Dn[j];
                        break;
                    case 2:
                        dop = d->Dn[k];
                        break;
                    }
                    d->Fmesh->fields[l][c_index]   = -2.0*PI*Leni*dop*d->Qmesh[c_index+1];
                    d->Fmesh->fields[l][c_index+1] =  2.0*PI*Leni*dop*d->Qmesh[c_index];
                }

            }

  #ifdef __detailed_timings
    timer = MPI_Wtime() - timer;
    t_convolution[1] = timer;
   timer = MPI_Wtime();
  #endif

    /* Durchfuehren der Fourier-Rueck-Transformation: */
    backward_fft(d);

  #ifdef __detailed_timings
    timer = MPI_Wtime() - timer;
    t_fft[1] += timer;
   timer = MPI_Wtime();
  #endif


    /* force assignment */
    assign_forces ( 1.0 / ( 2.0*s->length*s->length*s->length ), s, p, d, f, 0 );
    assign_forces ( 1.0 / ( 2.0*s->length*s->length*s->length ), s, p, d, f, 1 );

  #ifdef __detailed_timings
    timer = MPI_Wtime() - timer;
    t_force_assignment[1] = timer;
  #endif

    return;
}


void p3m_tune_aliasing_sums_ik_i (int nx, int ny, int nz, 
					   system_t *s, parameters_t *p, 
				  FLOAT_TYPE *alias1, FLOAT_TYPE *alias2, FLOAT_TYPE *alias3, FLOAT_TYPE *alias4)
{

  int    mx,my,mz;
  FLOAT_TYPE nmx,nmy,nmz;
  FLOAT_TYPE fnmx,fnmy,fnmz;

  FLOAT_TYPE ex,ex2,nm2,U2,factor1;
  int mesh = p->mesh;
  FLOAT_TYPE mesh_i = 1.0 / mesh;

  factor1 = SQR(PI / ( p->alpha * s->length ) );

  *alias1 = *alias2 = *alias3 = *alias4 = 0.0;
  for (mx=-P3M_BRILLOUIN_TUNING; mx<=P3M_BRILLOUIN_TUNING; mx++) {
    fnmx = mesh_i * (nmx = nx + mx*mesh);
    for (my=-P3M_BRILLOUIN_TUNING; my<=P3M_BRILLOUIN_TUNING; my++) {
      fnmy = mesh_i * (nmy = ny + my*mesh);
      for (mz=-P3M_BRILLOUIN_TUNING; mz<=P3M_BRILLOUIN_TUNING; mz++) {

	fnmz = mesh_i * (nmz = nz + mz*mesh);
	
	nm2 = SQR(nmx) + SQR(nmy) + SQR(nmz);
	ex = EXP(-factor1*nm2);
	ex2 = SQR( ex );
	
	U2 = my_power(sinc(fnmx)*sinc(fnmy)*sinc(fnmz), 2*p->cao);
	
	*alias1 += ex2 / nm2;
	*alias2 +=  U2 * ex * (nx*nmx + ny*nmy + nz*nmz) / nm2;
	*alias3 += U2;

       if (((mx+my+mz)%2)==0) {		//consider only even terms!
	 *alias4 += U2;
       } else {
	 *alias4 -= U2;
       }
      }
    }
  }
}



FLOAT_TYPE p3m_k_space_error_ik_i ( system_t *s, parameters_t *p ) {
    int  nx, ny, nz;
    FLOAT_TYPE he_q = 0.0;
    FLOAT_TYPE alias1, alias2, alias3, alias4, n2;
    int mesh = p->mesh;

    for ( nx=-mesh/2; nx<mesh/2; nx++ ) {
        for ( ny=-mesh/2; ny<mesh/2; ny++ ) {
            for ( nz=-mesh/2; nz<mesh/2; nz++ ) {
                if ( ( nx!=0 ) || ( ny!=0 ) || ( nz!=0 ) ) {
                    n2 = SQR ( nx ) + SQR ( ny ) + SQR ( nz );

                    p3m_tune_aliasing_sums_ik_i ( nx,ny,nz, s, p, &alias1,&alias2,&alias3, &alias4 );
                    he_q +=  alias1  -  SQR ( alias2 ) / (0.5*n2*(SQR(alias3)+SQR(alias4)) );
                }
            }
        }
    }
    he_q = FLOAT_ABS(he_q);
    return 2.0*s->q2*SQRT ( he_q/ ( FLOAT_TYPE ) s->nparticles ) / ( SQR ( s->length ) );
}

/*
FLOAT_TYPE p3m_k_space_error_ik_i ( system_t *s, parameters_t *p )
{
  int    NX,NY,NZ, N2;
  FLOAT_TYPE Dnx,Dny,Dnz;
  FLOAT_TYPE dMesh,dMeshi;
  FLOAT_TYPE Zaehler[3]={0.0,0.0,0.0},Nenner1=0.0, Nenner2=0.0;
  FLOAT_TYPE zwi;
  FLOAT_TYPE G_hat;
  FLOAT_TYPE alias1, alias2;
  FLOAT_TYPE he_q=0.0;
  FLOAT_TYPE fact =  p->mesh*p->mesh*p->mesh*2.0/(s->length*s->length);
  FLOAT_TYPE ctan_x, ctan_y, cs; 
  FLOAT_TYPE meshi = 1.0/(FLOAT_TYPE)(p->mesh);

  data_t d;
  
  int Mesh = p->mesh;
  FLOAT_TYPE Leni = 1.0/s->length;
  FLOAT_TYPE Aik = 1.0;

  d.nshift = Init_array( p->mesh, sizeof(FLOAT_TYPE) );
  d.Dn = Init_array ( p->mesh, sizeof(FLOAT_TYPE) );
  d.mesh = p->mesh;

  Init_nshift( &d );
  Init_differential_operator( &d );

  dMesh = (FLOAT_TYPE)Mesh;
  dMeshi= 1.0/dMesh;

  for (NX=0; NX<Mesh; NX++) {
    ctan_x = analytic_cotangent_sum ( d.nshift[NX], meshi,p->cao );
    for (NY=0; NY<Mesh; NY++) {
      ctan_y = ctan_x * analytic_cotangent_sum ( d.nshift[NY], meshi, p->cao );   
      for (NZ=0; NZ<Mesh; NZ++) {
	      
	if ((NX!=0) || (NY!=0) || (NZ!=0)) {
	  if ((NX%(Mesh/2) == 0) && (NY%(Mesh/2) == 0) && (NZ%(Mesh/2) == 0)) {
	    G_hat=0.0;
	    Aik = 1.0;
	  }
	  else {
	    Aliasing_sums_ik_i( s, p, &d, NX, NY, NZ, Zaehler, &Nenner1,  &Nenner2);

	    Dnx = d.Dn[NX];
	    Dny = d.Dn[NY];
	    Dnz = d.Dn[NZ];
	    
	    zwi  = Dnx*Zaehler[0] + Dny*Zaehler[1] + Dnz*Zaehler[2];
	    zwi /= ( SQR(Dnx) + SQR(Dny) + SQR(Dnz) );

	    Aik = 0.5*(SQR(Nenner1) + SQR(Nenner2));

	    zwi /= Aik;
	  
	    G_hat = zwi;
	  
	  }


	  p3m_tune_aliasing_sums_ik_i( d.nshift[NX], d.nshift[NY], d.nshift[NZ], s, p, &alias1, &alias2 );
	  
	  cs = ctan_y * analytic_cotangent_sum( d.nshift[NZ], meshi, p->cao );
	  
	  
	  N2 = SQR( d.nshift[NX] ) + SQR ( d.nshift[NY] ) + SQR ( d.nshift[NZ] );
		
	  // Neelov 2010, Eq. (B14)
	  //          he_q += (alias1  + SQR(cs)*N2*SQR(G_hat)/(fact*fact) -2.0*alias2*G_hat/fact);
          he_q += ( alias1  -  SQR ( alias2/cs ) / n2 );
	  //he_q += (alias1  + SQR(cs)*SQR(G_hat) - 2.0*alias2*G_hat );
	  //he_q += alias1 + SQR(G_hat) * ( Aik - N2*SQR(Leni*PI*2.0)*alias2 );
	  //	  printf("%d %d %d %e %e %e alias1 %e alias2 %e G_hat %e cs %e he_q %e fact %e\n", NX, NY, NZ, d.nshift[NX], d.nshift[NY], d.nshift[NZ], alias1, alias2, G_hat, cs, he_q, fact);
	}
      }
    }
  }
  FFTW_FREE( d.nshift );
  FFTW_FREE( d.Dn );

  he_q = fabs(he_q);                                                                                                                                                                                                                       
  return 2.0*s->q2*sqrt ( he_q/ ( FLOAT_TYPE ) s->nparticles ) / ( SQR ( s->length ) );  
}
*/


FLOAT_TYPE Error_ik_i( system_t *s, parameters_t *p) {
  FLOAT_TYPE real = Realspace_error( s, p);
  FLOAT_TYPE recp = p3m_k_space_error_ik_i( s, p );

  return SQRT( SQR ( real ) + SQR ( recp ) );
}

