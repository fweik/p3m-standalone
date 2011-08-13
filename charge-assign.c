#include "charge-assign.h"

#include <stdio.h>
#include <math.h>

// #define CA_DEBUG

#define r_ind(A,B,C) ((A)*Mesh*Mesh + (B)*Mesh + (C))
#define c_ind(A,B,C) (2*Mesh*Mesh*(A)+2*Mesh*(B)+2*(C))

void assign_charge(int id, FLOAT_TYPE q,
		       FLOAT_TYPE real_pos[3],
		   FLOAT_TYPE *p3m_rs_mesh, int ii)
{
  int d, i0, i1, i2;
  FLOAT_TYPE tmp0, tmp1;
  /* position of a particle in local mesh units */
  FLOAT_TYPE pos;
  /* 1d-index of nearest mesh point */
  int nmp;
  /* index for caf interpolation grid */
  int arg[3];
  /* index, index jumps for rs_mesh array */
  FLOAT_TYPE cur_ca_frac_val, *cur_ca_frac;
  int cf_cnt = id*cao3;

  int base[3];
  int i,j,k;
  FLOAT_TYPE MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;   

  FLOAT_TYPE Hi = (double)Mesh/(double)Len;

  int MESHMASKE = Mesh - 1;
  double pos_shift;

  pos_shift = (double)((cao-1)/2);
    /* particle position in mesh coordinates */
    for(d=0;d<3;d++) {
      pos    = real_pos[d]*Hi - pos_shift;
      nmp = (int) floor(pos + 0.5);
      base[d]  = (nmp > 0) ? nmp%MESHMASKE : (nmp + Mesh)%MESHMASKE;
      arg[d] = (int) floor((pos - nmp + 0.5)*MI2);
      ca_ind[ii][3*id + d] = base[d];
#ifdef CA_DEBUG
      fprintf(stderr, "%d: pos_shift %lf base_pos %lf, nmp %d, base %d, arg[d]*M2 %lf\n", d, pos_shift, real_pos[d] - pos_shift, nmp, base[d], arg[d]/MI2);
#endif
    }

    for(i0=0; i0<cao; i0++) {
      i = (base[0] + i0)&MESHMASKE;
      tmp0 = LadInt[i0][arg[0]];
      for(i1=0; i1<cao; i1++) {
        j = (base[1] + i1)&MESHMASKE;
	tmp1 = tmp0 * LadInt[i1][arg[1]];
	for(i2=0; i2<cao; i2++) {
          k = (base[2] + i2)&MESHMASKE;
	  cur_ca_frac_val = q * tmp1 * LadInt[i2][arg[2]];
          cf[ii][cf_cnt++] = cur_ca_frac_val ;
	  p3m_rs_mesh[c_ind(i,j,k)+ii] += cur_ca_frac_val;
	}
      }
    }
    //    printf("Assigned %lf total from particle %d\n", charge, id);
}

// assign the forces obtained from k-space 
void assign_forces(double force_prefac, FLOAT_TYPE *F, int Teilchenzahl, FLOAT_TYPE *p3m_rs_mesh,int ii) 
{
  int i,c,i0,i1,i2;
  int cf_cnt=0;
  int *base;
  int j,k,l;
  int MESHMASKE = Mesh - 1;
  FLOAT_TYPE A,B;

  cf_cnt=0;
  for(i=0; i<Teilchenzahl; i++) { 
    base = ca_ind[ii] + 3*i;
    for(i0=0; i0<cao; i0++) {
      j = (base[0] + i0)&MESHMASKE;
      for(i1=0; i1<cao; i1++) {
	k = (base[1] + i1)&MESHMASKE;
	for(i2=0; i2<cao; i2++) {
	  l = (base[2] + i2)&MESHMASKE;
	  A = cf[ii][cf_cnt];
	  B = p3m_rs_mesh[c_ind(j,k,l)+ii];

	  F[i] -= force_prefac*A*B; 
	  cf_cnt++;
	}
      }
    }
    if(ii==1)
      F[i] *= 0.5;
  }
}

void assign_charge_and_derivatives(int id, FLOAT_TYPE q,
		       FLOAT_TYPE real_pos[3],
				   FLOAT_TYPE *p3m_rs_mesh, int ii)
{
  int d, i0, i1, i2;
  FLOAT_TYPE tmp0, tmp1, tmp2;
  FLOAT_TYPE tmp0_x, tmp1_y, tmp2_z;
  /* position of a particle in local mesh units */
  FLOAT_TYPE pos;
  /* 1d-index of nearest mesh point */
  int nmp;
  /* index for caf interpolation grid */
  int arg[3];
  /* index, index jumps for rs_mesh array */
  FLOAT_TYPE cur_ca_frac_val, *cur_ca_frac;
  int cf_cnt = id*cao3;

  int base[3];
  int i,j,k;
  FLOAT_TYPE MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;   

  FLOAT_TYPE Hi = (double)Mesh/(double)Len;

  int MESHMASKE = Mesh - 1;
  FLOAT_TYPE charge = 0.0;
  double pos_shift;
  FLOAT_TYPE dMesh = (double)Mesh;

  pos_shift = (double)((cao-1)/2);
    /* particle position in mesh coordinates */
    for(d=0;d<3;d++) {
      pos    = real_pos[d]*Hi - pos_shift;
      nmp = (int) floor(pos + 0.5);
      base[d]  = (nmp > 0) ? nmp%MESHMASKE : (nmp + Mesh)%MESHMASKE;
      arg[d] = (int) floor((pos - nmp + 0.5)*MI2);
      ca_ind[ii][3*id + d] = base[d];
#ifdef CA_DEBUG
      fprintf(stderr, "%d: pos_shift %lf base_pos %lf, nmp %d, base %d, arg[d]*M2 %lf\n", d, pos_shift, real_pos[d] - pos_shift, nmp, base[d], arg[d]/MI2);
#endif
    }

    for(i0=0; i0<cao; i0++) {
      i = (base[0] + i0)&MESHMASKE;
      tmp0 = LadInt[i0][arg[0]];
      tmp0_x = LadInt_[i0][arg[0]]; 
      for(i1=0; i1<cao; i1++) {
        j = (base[1] + i1)&MESHMASKE;
	tmp1 = LadInt[i1][arg[1]];
        tmp1_y = LadInt_[i1][arg[1]];
	for(i2=0; i2<cao; i2++) {
          k = (base[2] + i2)&MESHMASKE;
          tmp2 = LadInt[i2][arg[2]];
          tmp2_z = LadInt_[i2][arg[2]];
	  cur_ca_frac_val = q * tmp0 * tmp1 * tmp2;
          cf[ii][cf_cnt] = cur_ca_frac_val ;
          dQdx[ii][cf_cnt] = tmp0_x * tmp1 * tmp2 * q;
          dQdy[ii][cf_cnt] = tmp0 * tmp1_y * tmp2 * q;
          dQdz[ii][cf_cnt] = tmp0 * tmp1 * tmp2_z * q;

	  p3m_rs_mesh[c_ind(i,j,k)+ii] += cur_ca_frac_val;
          cf_cnt++;
	}
      }
    }
}

// assign the forces obtained from k-space 
void assign_forces_ad(double force_prefac, FLOAT_TYPE **F, int Teilchenzahl, FLOAT_TYPE *p3m_rs_mesh,int ii) 
{
  int i,c,i0,i1,i2;
  int cf_cnt=0;
  int *base;
  int j,k,l;
  int MESHMASKE = Mesh - 1;
  FLOAT_TYPE A,B;

  cf_cnt=0;
  for(i=0; i<Teilchenzahl; i++) { 
    base = ca_ind[ii] + 3*i;
    cf_cnt = i*cao*cao*cao;
    for(i0=0; i0<cao; i0++) {
      j = (base[0] + i0)&MESHMASKE;
      for(i1=0; i1<cao; i1++) {
	k = (base[1] + i1)&MESHMASKE;
	for(i2=0; i2<cao; i2++) {
	  l = (base[2] + i2)&MESHMASKE;
	  B = p3m_rs_mesh[c_ind(j,k,l)+ii];
	  F[0][i] -= force_prefac*B*dQdx[ii][cf_cnt]; 
	  F[1][i] -= force_prefac*B*dQdy[ii][cf_cnt];
	  F[2][i] -= force_prefac*B*dQdz[ii][cf_cnt];
	  cf_cnt++;
	}
      }
    }
    if(ii==1) {
      F[0][i] *= 0.5;
      F[1][i] *= 0.5;
      F[2][i] *= 0.5;
    }
  }
}
