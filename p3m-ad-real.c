/**    Copyright (C) 2011,2012,2013 Florian Weik <fweik@icp.uni-stuttgart.de>

       This program is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

       This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this program.  If not, see <http://www.gnu.org/licenses/>. **/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <string.h>
#include <assert.h>

#include "types.h"
#include "charge-assign.h"
#include "p3m-ad-real.h"
#include "common.h"
#include "p3m-ad-self-forces.h"

#include "p3m-ad.h"

#include "realpart.h"

#ifdef __detailed_timings
#include <mpi.h>
#endif


const method_t method_p3m_ad_r = { METHOD_P3M_ad_r, "P3M with analytic differentiation, not intelaced, real input.", "p3m-ad-r",
				 METHOD_FLAG_P3M | METHOD_FLAG_ad | METHOD_FLAG_self_force_correction, 
				 &Init_ad_r, &Influence_function_berechnen_ad, &P3M_ad_r, &Error_ad, &p3m_k_space_error_ad };

static void forward_fft( data_t *d );
static void backward_fft( data_t *d );

inline void forward_fft ( data_t *d ) {
  FFTW_EXECUTE ( d->forward_plan[0] );
}

inline void backward_fft ( data_t *d ) {
  FFTW_EXECUTE ( d->backward_plan[0] );
}

data_t *Init_ad_r ( system_t *s, parameters_t *p ) {
    int mesh = p->mesh;

    data_t *d = Init_data ( &method_p3m_ad_r, s, p );

    d->forward_plans = 1;
    d->backward_plans = 1;

    d->forward_plan[0] = FFTW_PLAN_DFT_R2C_3D ( mesh, mesh, mesh, d->Qmesh, ( FFTW_COMPLEX * ) d->Qmesh, FFTW_PATIENT );

    d->backward_plan[0] = FFTW_PLAN_DFT_C2R_3D ( mesh, mesh, mesh, ( FFTW_COMPLEX * ) ( d->Qmesh ),  d->Qmesh,  FFTW_PATIENT );

    return d;
}

void P3M_ad_r( system_t *s, parameters_t *p, data_t *d, forces_t *f )
{
  
  /* Loop counters */
  int i, j, k, c_index; 
  /* Helper variables */
  FLOAT_TYPE T1;
  FLOAT_TYPE Leni = 1.0/s->length;
  int Mesh = p->mesh;
  
  memset(d->Qmesh, 0, 2*Mesh*Mesh*Mesh * sizeof(FLOAT_TYPE));

  TIMING_START_C
  
  /* chargeassignment */
  assign_charge_and_derivatives_real( s, p, d);

  TIMING_STOP_C
  TIMING_START_G
  
  /* Forward Fast Fourier Transform */
  forward_fft(d);

  for (i=0; i<Mesh; i++)
    for (j=0; j<Mesh; j++)
      for (k=0; k<(Mesh/2+1); k++)
	{
          c_index = 2*(Mesh*(Mesh/2+1) * i + (Mesh/2+1)*j + k);

	  T1 = d->G_hat[r_ind(i,j,k)];
	  d->Qmesh[c_index] *= T1;
	  d->Qmesh[c_index+1] *= T1;
	}

  /* Backward FFT */
  backward_fft(d);

  TIMING_STOP_G
  TIMING_START_F

  /* Force assignment */
  assign_forces_ad_real( Mesh * Leni * Leni * Leni , s, p, d, f);

#ifdef P3M_AD_SELF_FORCES
  Substract_self_forces(s,p,d,f);
#endif

  TIMING_STOP_F
  return;
}


