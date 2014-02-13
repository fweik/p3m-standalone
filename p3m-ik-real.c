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

// General typ definitions
#include "types.h"
#include "common.h"
#include "p3m-common.h"
// Charge assignment
#include "charge-assign.h"

// For realpart error
#include "realpart.h"

#include "p3m-ik-real.h"
#include "p3m-ik.h"

#ifdef __detailed_timings
#include <mpi.h>
#endif


// declaration of the method

const method_t method_p3m_ik_r = { METHOD_P3M_ik_r, "P3M with ik differentiation, not intelaced, real input.", "p3m-ik-r",
                                 METHOD_FLAG_P3M | METHOD_FLAG_ik,
                                 &Init_ik_r, &Influence_function_berechnen_ik, &P3M_ik_r, &Error_ik, &Error_ik_k,
                               };

// Forward declaration of local functions

static void forward_fft ( data_t * );
static void backward_fft ( data_t * );

inline void forward_fft ( data_t *d ) {
    FFTW_EXECUTE ( d->forward_plan[0] );
}

inline void backward_fft ( data_t *d ) {
    int i;
    for ( i=0;i<3;i++ )
        FFTW_EXECUTE ( d->backward_plan[i] );
}

data_t *Init_ik_r ( system_t *s, parameters_t *p ) {
    int l;
    int mesh = p->mesh;

    data_t *d = Init_data ( &method_p3m_ik_r, s, p );

    d->forward_plans = 1;
    d->backward_plans = 3;

    d->forward_plan[0] = FFTW_PLAN_DFT_R2C_3D ( mesh, mesh, mesh, d->Qmesh, (FFTW_COMPLEX *)d->Qmesh, FFTW_PATIENT );

    for ( l=0;l<3;l++ ) {
        d->backward_plan[l] = FFTW_PLAN_DFT_C2R_3D ( mesh, mesh, mesh, ( FFTW_COMPLEX * ) ( d->Fmesh->fields[l] ), ( d->Fmesh->fields[l] ), FFTW_PATIENT );
    }
    return d;
}

/* Calculates k-space part of the force, using ik-differentiation.
 */

void P3M_ik_r ( system_t *s, parameters_t *p, data_t *d, forces_t *f ) {
    /* Loop counters */
    int i, j, k;
    /* helper variables */
    FLOAT_TYPE T1;
    FLOAT_TYPE dop;

    // One over boxlength
    FLOAT_TYPE Leni = 1.0/s->length;

    int Mesh = p->mesh;
    int c_index;

    /* Setting charge mesh to zero */
    memset ( d->Qmesh, 0, 2*Mesh*Mesh*Mesh*sizeof ( FLOAT_TYPE ) );

    TIMING_START_C

    /* chargeassignment */
    assign_charge_real ( s, p, d );

    TIMING_STOP_C
    TIMING_START_G

    /* Forward Fast Fourier Transform */
    forward_fft(d);

    double q_r, q_i;

    /* Convolution */
    for ( i=0; i<Mesh; i++ ) {
      for ( j=0; j<Mesh; j++ ) {
	for ( k=0; k<(Mesh/2+1); k++ ) {
	  c_index = 2*(Mesh*(Mesh/2+1)* i + (Mesh/2+1)*j + k);

	  T1 = d->G_hat[r_ind ( i,j,k ) ];

	  q_r = d->Qmesh[c_index] *T1;
	  q_i = d->Qmesh[c_index+1] *T1;

	  dop = d->Dn[i];	 
	  d->Fmesh->fields[0][c_index]   =  -2.0*PI*Leni*dop*q_i;
	  d->Fmesh->fields[0][c_index+1] =   2.0*PI*Leni*dop*q_r;

	  dop = d->Dn[j];
	  d->Fmesh->fields[1][c_index]   =  -2.0*PI*Leni*dop*q_i;
	  d->Fmesh->fields[1][c_index+1] =   2.0*PI*Leni*dop*q_r;

	  dop = d->Dn[k];
	  d->Fmesh->fields[2][c_index]   =  -2.0*PI*Leni*dop*q_i;
	  d->Fmesh->fields[2][c_index+1] =   2.0*PI*Leni*dop*q_r;
	}
      }
    }
    /* Backward Fast Fourier Transformation */
    backward_fft(d);

    TIMING_STOP_G
    TIMING_START_F

    /* Force assignment */
    assign_forces_real ( 1.0/ ( 2.0*s->length*s->length*s->length ),s,p,d,f);

    TIMING_STOP_F
}

