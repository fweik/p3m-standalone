#include "charge-assign.h"

void assign_charge(int id, FLOAT_TYPE q,
		       FLOAT_TYPE real_pos[3],
		       FLOAT_TYPE *p3m_rs_mesh)
{
  int cao = ip+1;

  int d, i0, i1, i2;
  FLOAT_TYPE tmp0, tmp1;
  /* position of a particle in local mesh units */
  FLOAT_TYPE pos;
  /* 1d-index of nearest mesh point */
  int nmp;
  /* distance to nearest mesh point */
  FLOAT_TYPE dist[3];
  /* index for caf interpolation grid */
  int arg[3];
  /* index, index jumps for rs_mesh array */
  int q_ind = 0;
  FLOAT_TYPE cur_ca_frac_val, *cur_ca_frac;
  int cf_cnt = 0;

  FLOAT_TYPE MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;   

  FLOAT_TYPE Hi = (double)Mesh/(double)Len;

  int q_m_off = (Mesh - cao);
  int q_s_off = Mesh * q_m_off;

    /* particle position in mesh coordinates */
    for(d=0;d<3;d++) {
      pos    = real_pos[d]*Hi;
      nmp    = (int) pos;
      arg[d] = (int) ((pos - nmp)*MI2);
      /* for the first dimension, q_ind is always zero, so this shifts correctly */
      q_ind = nmp + Mesh*q_ind;
    }
    ca_ind[id] = q_ind;

    for(i0=0; i0<cao; i0++) {
      tmp0 = LadInt[i0][arg[0]];
      for(i1=0; i1<cao; i1++) {
	tmp1 = tmp0 * LadInt[i1][arg[1]];
	for(i2=0; i2<cao; i2++) {
	  cur_ca_frac_val = q * tmp1 * LadInt[i2][arg[2]];
          cf[cf_cnt++] = cur_ca_frac_val;
	  p3m_rs_mesh[2*q_ind] += cur_ca_frac_val;
	  q_ind++;
	}
	q_ind += q_m_off;
      }
      q_ind += q_s_off;
    }
  
}

// assign the forces obtained from k-space 
void assign_forces(double force_prefac, FLOAT_TYPE *F, int Teilchenzahl, FLOAT_TYPE *p3m_rs_mesh) 
{
  int i,c,i0,i1,i2;
  int cao = ip  + 1;
  int cp_cnt=0, cf_cnt=0;

  int q_ind;
  int q_m_off = (Mesh - cao);
  int q_s_off = Mesh * q_m_off;

  double A,B,C;
 
  cf_cnt=0;
    for(i=0; i<Teilchenzahl; i++) { 
      q_ind = ca_ind[i];
      for(i0=0; i0<cao; i0++) {
        for(i1=0; i1<cao; i1++) {
          for(i2=0; i2<cao; i2++) {
            A = cf[cf_cnt];
            B = p3m_rs_mesh[2*q_ind++];
            F[i] -= force_prefac*A*B; 
	    cf_cnt++;
	  }
	  q_ind += q_m_off;
	}
	q_ind += q_s_off;
      }
    }
}


