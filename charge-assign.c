void p3m_assign_charge(double q,
		       double real_pos[3],
		       int cp_cnt)
{
  extern double p3m_caf(int i, double xc,int cao_value);
  extern void p3m_realloc_ca_fields(int size);

  extern int    *ca_fmp;
  extern double *ca_frac;
  extern double *int_caf[7];
  extern double *p3m_rs_mesh;

  int d, i0, i1, i2;
  double tmp0, tmp1;
  /* position of a particle in local mesh units */
  double pos;
  /* 1d-index of nearest mesh point */
  int nmp;
  /* distance to nearest mesh point */
  double dist[3];
  /* index for caf interpolation grid */
  int arg[3];
  /* index, index jumps for rs_mesh array */
  int q_ind = 0;
  double cur_ca_frac_val, *cur_ca_frac;

  // make sure we have enough space
  if (cp_cnt >= ca_num) p3m_realloc_ca_fields(cp_cnt + 1);
  // do it here, since p3m_realloc_ca_fields may change the address of ca_frac
  cur_ca_frac = ca_frac + p3m.cao3*cp_cnt;

  if (p3m.inter == 0) {
    for(d=0;d<3;d++) {
      /* particle position in mesh coordinates */
      pos    = ((real_pos[d]-p3m_lm.ld_pos[d])*p3m.ai[d]);
      /* nearest mesh point */
      nmp  = (int)pos;
      /* distance to nearest mesh point */
      dist[d] = (pos-nmp)-0.5;
      /* 3d-array index of nearest mesh point */
      q_ind = (d == 0) ? nmp : nmp + p3m_lm.dim[d]*q_ind;
    }
    if (cp_cnt >= 0) ca_fmp[cp_cnt] = q_ind;
    
    for(i0=0; i0<p3m.cao; i0++) {
      tmp0 = p3m_caf(i0, dist[0],p3m.cao);
      for(i1=0; i1<p3m.cao; i1++) {
	tmp1 = tmp0 * p3m_caf(i1, dist[1],p3m.cao);
	for(i2=0; i2<p3m.cao; i2++) {
	  cur_ca_frac_val = q * tmp1 * p3m_caf(i2, dist[2], p3m.cao);
	  if (cp_cnt >= 0) *(cur_ca_frac++) = cur_ca_frac_val;
	  p3m_rs_mesh[q_ind] += cur_ca_frac_val;
	  q_ind++;
	}
	q_ind += p3m_lm.q_2_off;
      }
      q_ind += p3m_lm.q_21_off;
    }
  }
  else {
    /* particle position in mesh coordinates */
    for(d=0;d<3;d++) {
      pos    = ((real_pos[d]-p3m_lm.ld_pos[d])*p3m.ai[d]);
      nmp    = (int) pos;
      arg[d] = (int) ((pos - nmp)*p3m.inter2);
      /* for the first dimension, q_ind is always zero, so this shifts correctly */
      q_ind = nmp + p3m_lm.dim[d]*q_ind;
    }
    if (cp_cnt >= 0) ca_fmp[cp_cnt] = q_ind;

    for(i0=0; i0<p3m.cao; i0++) {
      tmp0 = int_caf[i0][arg[0]];
      for(i1=0; i1<p3m.cao; i1++) {
	tmp1 = tmp0 * int_caf[i1][arg[1]];
	for(i2=0; i2<p3m.cao; i2++) {
	  cur_ca_frac_val = q * tmp1 * int_caf[i2][arg[2]];
	  if (cp_cnt >= 0) *(cur_ca_frac++) = cur_ca_frac_val;
	  p3m_rs_mesh[q_ind] += cur_ca_frac_val;
	  q_ind++;
	}
	q_ind += p3m_lm.q_2_off;
      }
      q_ind += p3m_lm.q_21_off;
    }
  }
}

