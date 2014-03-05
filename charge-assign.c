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

#include "charge-assign.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>

// #define CA_DEBUG

inline static int wrap_mesh_index(int ind, int mesh) {
  if((ind > 0) && (ind < mesh))
    return ind;
  else if(ind < 0)
    ind += mesh;
  else if(ind >= mesh)
    ind -= mesh;

  // mesh > cao/2 !
  /* if((ret < 0) || (ret >= mesh)) */
  /*   return wrap_mesh_index(ret, mesh); */
  return ind;
}

inline static int int_floor(FLOAT_TYPE x)
{
  int i = (int)x; /* truncate */
  int n = ( x != (FLOAT_TYPE)i );
  int g = ( x < 0 );
  return i - ( n & g ); /* i-1 if x<0 and x!=i */
}

void assign_charge(system_t *s, parameters_t *p, data_t *d, int ii)
{
    int dim, i0, i1, i2, id;
    FLOAT_TYPE tmp0, tmp1;
    /* position of a particle in local mesh units */
    FLOAT_TYPE pos;
    /* 1d-index of nearest mesh point */
    int nmp;
    /* index for caf interpolation grid */
    int arg[3];
    /* index, index jumps for rs_mesh array */
    FLOAT_TYPE cur_ca_frac_val;
    FLOAT_TYPE *cf_cnt;
    // Mesh coordinates of the closest mesh point
    int base[3];
    int i,j,k;
    FLOAT_TYPE MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;

    FLOAT_TYPE Hi = (double)d->mesh/(double)s->length;

    FLOAT_TYPE *cf = d->cf[ii];
    FLOAT_TYPE **interpol = d->inter->interpol;
    FLOAT_TYPE *Qmesh = d->Qmesh;
    FLOAT_TYPE q;
    const int cao = p->cao;
    const int mesh = d->mesh;

    // Make sure parameter-set and data-set are compatible

    FLOAT_TYPE pos_shift;

    /* Shift for odd charge assignment order */
    pos_shift = (FLOAT_TYPE)((p->cao-1)/2);

    for (id=0;id<s->nparticles;id++) {
        /* particle position in mesh coordinates */
        for (dim=0;dim<3;dim++) {
            pos    = s->p->fields[dim][id]*Hi - pos_shift + 0.5*ii;
            nmp = int_floor(pos + 0.5);
	    base[dim]  = wrap_mesh_index( nmp, d->mesh);
            arg[dim] = int_floor((pos - nmp + 0.5)*MI2);
            d->ca_ind[ii][3*id + dim] = base[dim];
        }
	q = s->q[id];
        cf_cnt = cf + id*p->cao3;
	for (i0=0; i0<cao; i0++) {
	  i = wrap_mesh_index(base[0] + i0, mesh);
	  tmp0 = q * interpol[arg[0]][i0];
	  for (i1=0; i1<cao; i1++) {
	    tmp1 = tmp0 * interpol[arg[1]][i1];
	    j = wrap_mesh_index(base[1] + i1, mesh);
	    for (i2=0; i2<cao; i2++) {
	      cur_ca_frac_val = tmp1 * interpol[arg[2]][i2];
	      k = wrap_mesh_index(base[2] + i2, mesh);
	      *cf_cnt++ = cur_ca_frac_val;
	      Qmesh[c_ind(i,j,k)+ii] += cur_ca_frac_val;
	    }
	  }
	} 
    }
}

void assign_charge_real(system_t *s, parameters_t *p, data_t *d)
{
    int dim, i0, i1, i2, id;
    FLOAT_TYPE tmp0, tmp1;
    /* position of a particle in local mesh units */
    FLOAT_TYPE pos;
    /* 1d-index of nearest mesh point */
    int nmp;
    /* index for caf interpolation grid */
    int arg[3];
    /* index, index jumps for rs_mesh array */
    FLOAT_TYPE cur_ca_frac_val;
    FLOAT_TYPE *cf_cnt;
    // Mesh coordinates of the closest mesh point
    int base[3];
    int i,j,k;
    FLOAT_TYPE MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;

    FLOAT_TYPE Hi = (double)d->mesh/(double)s->length;

    FLOAT_TYPE *cf = d->cf[0];
    FLOAT_TYPE **interpol = d->inter->interpol;
    FLOAT_TYPE *Qmesh = d->Qmesh;
    FLOAT_TYPE q;
    const int cao = p->cao;
    const int mesh = d->mesh;

    // Make sure parameter-set and data-set are compatible

    FLOAT_TYPE pos_shift;

    /* Shift for odd charge assignment order */
    pos_shift = (FLOAT_TYPE)((p->cao-1)/2);

    for (id=0;id<s->nparticles;id++) {
        /* particle position in mesh coordinates */
        for (dim=0;dim<3;dim++) {
            pos    = s->p->fields[dim][id]*Hi - pos_shift;
            nmp = int_floor(pos + 0.5);
	    base[dim]  = wrap_mesh_index( nmp, d->mesh);
            arg[dim] = int_floor((pos - nmp + 0.5)*MI2);
            d->ca_ind[0][3*id + dim] = base[dim];
        }
	q = s->q[id];
        cf_cnt = cf + id*p->cao3;
	for (i0=0; i0<cao; i0++) {
	  i = wrap_mesh_index(base[0] + i0, mesh);
	  tmp0 = q * interpol[arg[0]][i0];
	  for (i1=0; i1<cao; i1++) {
	    tmp1 = tmp0 * interpol[arg[1]][i1];
	    j = wrap_mesh_index(base[1] + i1, mesh);
	    for (i2=0; i2<cao; i2++) {
	      cur_ca_frac_val = tmp1 * interpol[arg[2]][i2];
	      k = wrap_mesh_index(base[2] + i2, mesh);
	      *cf_cnt++ = cur_ca_frac_val;
	      Qmesh[mesh*(mesh+2) * i + (mesh+2) * j + k] += cur_ca_frac_val;
	    }
	  }
	} 
    }
}


void assign_charge_q2(system_t *s, parameters_t *p, FLOAT_TYPE *Qmesh, int mesh, interpolation_t *inter)
{
    int dim, i0, i1, i2, id;
    FLOAT_TYPE tmp0, tmp1;
    /* position of a particle in local mesh units */
    FLOAT_TYPE pos;
    /* 1d-index of nearest mesh point */
    int nmp;
    /* index for caf interpolation grid */
    int arg[3];
    /* index, index jumps for rs_mesh array */
    FLOAT_TYPE cur_ca_frac_val;
    // Mesh coordinates of the closest mesh point
    int base[3];
    int i,j,k;
    FLOAT_TYPE MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;

    FLOAT_TYPE Hi = (double)p->mesh/(double)s->length;

    FLOAT_TYPE **interpol = inter->interpol;
    FLOAT_TYPE q2;
    const int cao = p->cao;

    // Make sure parameter-set and data-set are compatible

    FLOAT_TYPE pos_shift;

    /* Shift for odd charge assignment order */
    pos_shift = (FLOAT_TYPE)((p->cao-1)/2);

    for (id=0;id<s->nparticles;id++) {
        /* particle position in mesh coordinates */
        for (dim=0;dim<3;dim++) {
            pos    = s->p->fields[dim][id]*Hi - pos_shift;
            nmp = int_floor(pos + 0.5);
	    base[dim]  = wrap_mesh_index( nmp, mesh);
            arg[dim] = int_floor((pos - nmp + 0.5)*MI2);
        }
	q2 = SQR(s->q[id]);
	for (i0=0; i0<cao; i0++) {
	  i = wrap_mesh_index(base[0] + i0, mesh);
	  tmp0 = q2 * interpol[arg[0]][i0];
	  for (i1=0; i1<cao; i1++) {
	    tmp1 = tmp0 * interpol[arg[1]][i1];
	    j = wrap_mesh_index(base[1] + i1, mesh);
	    for (i2=0; i2<cao; i2++) {
	      cur_ca_frac_val = tmp1 * interpol[arg[2]][i2];
	      k = wrap_mesh_index(base[2] + i2, mesh);
	      Qmesh[2*(mesh*mesh*i+mesh*j+k)] += cur_ca_frac_val;
	    }
	  }
	} 
    }
}

/* // assign the forces obtained from k-space */
/* void assign_q2(FLOAT_TYPE force_prefac, system_t *s, parameters_t *p, FLOAT_TYPE *Fmesh, forces_t *f) { */
/*   int i,i0,i1,i2; */
/*   FLOAT_TYPE *cf_cnt; */
/*   int *base; */
/*   int j,k,l; */
/*   FLOAT_TYPE B; */
/*   FLOAT_TYPE field_x, field_y, field_z; */
/*   int l_ind; */
/*   FLOAT_TYPE *fmesh_x = d->Fmesh->fields[0], *fmesh_y = d->Fmesh->fields[1], *fmesh_z = d->Fmesh->fields[2]; */
/*   int mesh = d->mesh; */

/*   const int cao = p->cao; */

/*   for (id=0;id<s->nparticles;id++) { */
/*     /\* particle position in mesh coordinates *\/ */
/*     for (dim=0;dim<3;dim++) { */
/*       pos    = s->p->fields[dim][id]*Hi - pos_shift; */
/*       nmp = int_floor(pos + 0.5); */
/*       base[dim]  = wrap_mesh_index( nmp, mesh); */
/*       arg[dim] = int_floor((pos - nmp + 0.5)*MI2); */
/*     } */
/*     q2 = SQR(s->q[id]); */

/*   for (i=0; i<s->nparticles; i++) { */
/*     field_x = field_y = field_z = 0; */
/*     base = d->ca_ind[ii] + 3*i; */
/*     for (i0=0; i0<cao; i0++) { */
/*       j = wrap_mesh_index(base[0] + i0, mesh); */
/*       for (i1=0; i1<cao; i1++) { */
/* 	k = wrap_mesh_index(base[1] + i1, mesh); */
/* 	for (i2=0; i2<cao; i2++) { */
/* 	  l = wrap_mesh_index(base[2] + i2, mesh); */
/* 	  B = force_prefac*(*cf_cnt++); */
/* 	  l_ind = c_ind(j,k,l)+ii; */
/* 	  field_x -= fmesh_x[l_ind]*B; */
/* 	  field_y -= fmesh_y[l_ind]*B; */
/* 	  field_z -= fmesh_z[l_ind]*B; */
/* 	} */
/*       } */
/*     } */
/*     f->f_k->fields[0][i] += field_x; */
/*     f->f_k->fields[1][i] += field_y; */
/*     f->f_k->fields[2][i] += field_z; */
/*   } */

/* } */


void assign_charge_nocf(system_t *s, parameters_t *p, FLOAT_TYPE *Qmesh, int mesh, interpolation_t *inter)
{
  int dim, i0, i1, i2, id;
  FLOAT_TYPE tmp0, tmp1;
  /* position of a particle in local mesh units */
  FLOAT_TYPE pos;
  /* 1d-index of nearest mesh point */
  int nmp;
  /* index for caf interpolation grid */
  int arg[3];
  /* index, index jumps for rs_mesh array */
  FLOAT_TYPE cur_ca_frac_val;
  // Mesh coordinates of the closest mesh point
  int base[3];
  int i,j,k;
  FLOAT_TYPE MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;

  FLOAT_TYPE Hi = (double)p->mesh/(double)s->length;

  FLOAT_TYPE **interpol = inter->interpol;
  FLOAT_TYPE q;
  const int cao = p->cao;

  // Make sure parameter-set and data-set are compatible

  FLOAT_TYPE pos_shift;

  /* Shift for odd charge assignment order */
  pos_shift = (FLOAT_TYPE)((p->cao-1)/2);

  for (id=0;id<s->nparticles;id++) {
    /* particle position in mesh coordinates */
    for (dim=0;dim<3;dim++) {
      pos    = s->p->fields[dim][id]*Hi - pos_shift;
      nmp = int_floor(pos + 0.5);
      base[dim]  = wrap_mesh_index( nmp, mesh);
      arg[dim] = int_floor((pos - nmp + 0.5)*MI2);
    }
    q = s->q[id];
    for (i0=0; i0<cao; i0++) {
      i = wrap_mesh_index(base[0] + i0, mesh);
      tmp0 = q * interpol[arg[0]][i0];
      for (i1=0; i1<cao; i1++) {
	tmp1 = tmp0 * interpol[arg[1]][i1];
	j = wrap_mesh_index(base[1] + i1, mesh);
	for (i2=0; i2<cao; i2++) {
	  cur_ca_frac_val = tmp1 * interpol[arg[2]][i2];
	  k = wrap_mesh_index(base[2] + i2, mesh);
	  /* printf("assignment %d %d %d cf %e\n", i, j, k, cur_ca_frac_val); */
	  Qmesh[2*(mesh*mesh*i+mesh*j+k)] += cur_ca_frac_val;
	}
      }
    } 
  }
}

void collect_rms_nocf(system_t *s, parameters_t *p, FLOAT_TYPE *Qmesh, FLOAT_TYPE *rms, int mesh, interpolation_t *inter)
{
  int dim, i0, i1, i2, id;
  FLOAT_TYPE tmp0, tmp1;
  /* position of a particle in local mesh units */
  FLOAT_TYPE pos;
  /* 1d-index of nearest mesh point */
  int nmp;
  /* index for caf interpolation grid */
  int arg[3];
  /* index, index jumps for rs_mesh array */
  FLOAT_TYPE cur_ca_frac_val;
  // Mesh coordinates of the closest mesh point
  int base[3];
  int i,j,k;
  FLOAT_TYPE MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;

  FLOAT_TYPE Hi = (double)p->mesh/(double)s->length;

  FLOAT_TYPE **interpol = inter->interpol;
  FLOAT_TYPE q2;
  const int cao = p->cao;

  // Make sure parameter-set and data-set are compatible

  FLOAT_TYPE pos_shift;

  /* Shift for odd charge assignment order */
  pos_shift = (FLOAT_TYPE)((p->cao-1)/2);

  for (id=0;id<s->nparticles;id++) {
    /* particle position in mesh coordinates */
    for (dim=0;dim<3;dim++) {
      pos    = s->p->fields[dim][id]*Hi - pos_shift;
      nmp = int_floor(pos + 0.5);
      base[dim]  = wrap_mesh_index( nmp, mesh);
      arg[dim] = int_floor((pos - nmp + 0.5)*MI2);
    }
    q2 = SQR(s->q[id]);
    rms[id] = 0.0;
    for (i0=0; i0<cao; i0++) {
      i = wrap_mesh_index(base[0] + i0, mesh);
      tmp0 = q2 * interpol[arg[0]][i0];
      for (i1=0; i1<cao; i1++) {
	tmp1 = tmp0 * interpol[arg[1]][i1];
	j = wrap_mesh_index(base[1] + i1, mesh);
	for (i2=0; i2<cao; i2++) {
	  cur_ca_frac_val = tmp1 * interpol[arg[2]][i2];
	  k = wrap_mesh_index(base[2] + i2, mesh);
	  /* printf("assignment %d %d %d cf %e\n", i, j, k, cur_ca_frac_val); */
	  rms[id] += Qmesh[2*(mesh*mesh*i+mesh*j+k)] * cur_ca_frac_val;
	}
      }
    } 
  }
}

// assign the forces obtained from k-space
void assign_forces(FLOAT_TYPE force_prefac, system_t *s, parameters_t *p, data_t *d, forces_t *f, int ii) {
  int i,i0,i1,i2;
  FLOAT_TYPE *cf_cnt;
  int *base;
  int j,k,l;
  FLOAT_TYPE B;
  FLOAT_TYPE field_x, field_y, field_z;
  int l_ind;
  FLOAT_TYPE *fmesh_x = d->Fmesh->fields[0], *fmesh_y = d->Fmesh->fields[1], *fmesh_z = d->Fmesh->fields[2];
  int mesh = d->mesh;

  const int cao = p->cao;

  cf_cnt = d->cf[ii];

  for (i=0; i<s->nparticles; i++) {
    field_x = field_y = field_z = 0;
    base = d->ca_ind[ii] + 3*i;
    for (i0=0; i0<cao; i0++) {
      j = wrap_mesh_index(base[0] + i0, mesh);
      for (i1=0; i1<cao; i1++) {
	k = wrap_mesh_index(base[1] + i1, mesh);
	for (i2=0; i2<cao; i2++) {
	  l = wrap_mesh_index(base[2] + i2, mesh);
	  B = force_prefac*(*cf_cnt++);
	  l_ind = c_ind(j,k,l)+ii;
	  field_x -= fmesh_x[l_ind]*B;
	  field_y -= fmesh_y[l_ind]*B;
	  field_z -= fmesh_z[l_ind]*B;
	}
      }
    }
    f->f_k->fields[0][i] += field_x;
    f->f_k->fields[1][i] += field_y;
    f->f_k->fields[2][i] += field_z;

    if (ii==1) {
      f->f_k->fields[0][i] *= 0.5;
      f->f_k->fields[1][i] *= 0.5;
      f->f_k->fields[2][i] *= 0.5;
    }
  }

}

// assign the forces obtained from k-space
void assign_forces_interlacing(FLOAT_TYPE force_prefac, system_t *s, parameters_t *p, data_t *d, forces_t *f) {
  int i,i0,i1,i2;
  FLOAT_TYPE *cf_cnt1, *cf_cnt2;
  int *base1, *base2;
  int j1,k1,l1,j2,k2,l2;
  FLOAT_TYPE B1, B2;
  FLOAT_TYPE field_x, field_y, field_z;
  int l_ind1, l_ind2;
  FLOAT_TYPE *fmesh_x = d->Fmesh->fields[0], *fmesh_y = d->Fmesh->fields[1], *fmesh_z = d->Fmesh->fields[2];
  int mesh = d->mesh;

  const int cao = p->cao;

  cf_cnt1 = d->cf[0];
  cf_cnt2 = d->cf[1];

  for (i=0; i<s->nparticles; i++) {
    field_x = field_y = field_z = 0;
    base1 = d->ca_ind[0] + 3*i;
    base2 = d->ca_ind[1] + 3*i;
    for (i0=0; i0<cao; i0++) {
      j1 = wrap_mesh_index(base1[0] + i0, mesh);
      j2 = wrap_mesh_index(base2[0] + i0, mesh);
      for (i1=0; i1<cao; i1++) {
	k1 = wrap_mesh_index(base1[1] + i1, mesh);
	k2 = wrap_mesh_index(base2[1] + i1, mesh);
	for (i2=0; i2<cao; i2++) {
	  l1 = wrap_mesh_index(base1[2] + i2, mesh);
	  l2 = wrap_mesh_index(base2[2] + i2, mesh);
	  B1 = 0.5*force_prefac*(*cf_cnt1++);
	  B2 = 0.5*force_prefac*(*cf_cnt2++);
	  l_ind1 = c_ind(j1,k1,l1);
	  l_ind2 = c_ind(j2,k2,l2);
	  field_x -= fmesh_x[l_ind1]*B1 + fmesh_x[l_ind2+1]*B2;
	  field_y -= fmesh_y[l_ind1]*B1 + fmesh_y[l_ind2+1]*B2;
	  field_z -= fmesh_z[l_ind1]*B1 + fmesh_z[l_ind2+1]*B2;
	}
      }
    }
    f->f_k->fields[0][i] += field_x;
    f->f_k->fields[1][i] += field_y;
    f->f_k->fields[2][i] += field_z;

  }

}

void assign_forces_interlacing_ad(FLOAT_TYPE force_prefac, system_t *s, parameters_t *p, data_t *d, forces_t *f) {
  int i,i0,i1,i2;
  int cf_cnt;
  int *base1, *base2;
  int j1,k1,l1,j2,k2,l2;
  FLOAT_TYPE B1, B2;
  FLOAT_TYPE field_x, field_y, field_z;
  int mesh = d->mesh;

  FLOAT_TYPE *dQdx0 = d->dQdx[0];
  FLOAT_TYPE *dQdy0 = d->dQdy[0];
  FLOAT_TYPE *dQdz0 = d->dQdz[0];

  FLOAT_TYPE *dQdx1 = d->dQdx[1];
  FLOAT_TYPE *dQdy1 = d->dQdy[1];
  FLOAT_TYPE *dQdz1 = d->dQdz[1];

  const int cao = p->cao;

  cf_cnt = 0;

  for (i=0; i<s->nparticles; i++) {
    field_x = field_y = field_z = 0;
    base1 = d->ca_ind[0] + 3*i;
    base2 = d->ca_ind[1] + 3*i;
    cf_cnt = i*p->cao3;
    for (i0=0; i0<cao; i0++) {
      j1 = wrap_mesh_index(base1[0] + i0, mesh);
      j2 = wrap_mesh_index(base2[0] + i0, mesh);
      for (i1=0; i1<cao; i1++) {
	k1 = wrap_mesh_index(base1[1] + i1, mesh);
	k2 = wrap_mesh_index(base2[1] + i1, mesh);
	for (i2=0; i2<cao; i2++) {
	  l1 = wrap_mesh_index(base1[2] + i2, mesh);
	  l2 = wrap_mesh_index(base2[2] + i2, mesh);
	  B1 = d->Qmesh[c_ind(j1,k1,l1)+0];
	  B2 = d->Qmesh[c_ind(j2,k2,l2)+1];
	  field_x -= 0.5*force_prefac*(B1*dQdx0[cf_cnt] + B2*dQdx1[cf_cnt]);
	  field_y -= 0.5*force_prefac*(B1*dQdy0[cf_cnt] + B2*dQdy1[cf_cnt]);
	  field_z -= 0.5*force_prefac*(B1*dQdz0[cf_cnt] + B2*dQdz1[cf_cnt]);
	  cf_cnt++;
	}
      }
    }
    f->f_k->fields[0][i] += field_x;
    f->f_k->fields[1][i] += field_y;
    f->f_k->fields[2][i] += field_z;
  }

}




void assign_forces_real(FLOAT_TYPE force_prefac, system_t *s, parameters_t *p, data_t *d, forces_t *f) {
  int i,i0,i1,i2;
  FLOAT_TYPE *cf_cnt;
  int *base;
  int j,k,l;
  FLOAT_TYPE B;
  FLOAT_TYPE field_x, field_y, field_z;
  int l_ind;
  FLOAT_TYPE *fmesh_x = d->Fmesh->fields[0], *fmesh_y = d->Fmesh->fields[1], *fmesh_z = d->Fmesh->fields[2];
  int mesh = d->mesh;

  const int cao = p->cao;

  cf_cnt = d->cf[0];

  for (i=0; i<s->nparticles; i++) {
    field_x = field_y = field_z = 0;
    base = d->ca_ind[0] + 3*i;
    for (i0=0; i0<cao; i0++) {
      j = wrap_mesh_index(base[0] + i0, mesh);
      for (i1=0; i1<cao; i1++) {
	k = wrap_mesh_index(base[1] + i1, mesh);
	for (i2=0; i2<cao; i2++) {
	  l = wrap_mesh_index(base[2] + i2, mesh);
	  B = force_prefac*(*cf_cnt++);
	  l_ind = mesh*(mesh+2) * j + (mesh+2) * k + l;
	  field_x -= fmesh_x[l_ind]*B;
	  field_y -= fmesh_y[l_ind]*B;
	  field_z -= fmesh_z[l_ind]*B;
	}
      }
    }
    f->f_k->fields[0][i] += field_x;
    f->f_k->fields[1][i] += field_y;
    f->f_k->fields[2][i] += field_z;

  }

}


void assign_charge_and_derivatives(system_t *s, parameters_t *p, data_t *d, int ii)
{
    int dim, i0, i1, i2;
    int id;
    FLOAT_TYPE tmp0, tmp1, tmp2;
    FLOAT_TYPE tmp0_x, tmp1_y, tmp2_z;
    /* position of a particle in local mesh units */
    FLOAT_TYPE pos;
    /* 1d-index of nearest mesh point */
    int nmp;
    /* index for caf interpolation grid */
    int arg[3];
    /* index, index jumps for rs_mesh array */
    FLOAT_TYPE cur_ca_frac_val;
    int cf_cnt;

    int base[3];
    int i,j,k;
    int Mesh = p->mesh;
    FLOAT_TYPE MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;

    FLOAT_TYPE Hi = (double)Mesh/s->length;
    FLOAT_TYPE Leni = 1.0/s->length;
    FLOAT_TYPE q, qLeni;

    double pos_shift;

    pos_shift = (double)((p->cao-1)/2);
    /* particle position in mesh coordinates */
    for (id=0;id<s->nparticles;id++) {
        cf_cnt = id*p->cao3;
	q = s->q[id] ;
	qLeni = q * Leni;
        for (dim=0;dim<3;dim++) {
	  pos    = s->p->fields[dim][id]*Hi - pos_shift + 0.5*ii;
	  nmp = int_floor(pos + 0.5);
	  base[dim]  = wrap_mesh_index( nmp, d->mesh);
	  arg[dim] = int_floor((pos - nmp + 0.5)*MI2);
	  d->ca_ind[ii][3*id + dim] = base[dim];
        }

        for (i0=0; i0<p->cao; i0++) {
	  i = wrap_mesh_index(base[0] + i0, d->mesh);
	  tmp0 = d->inter->interpol[arg[0]][i0];
	  tmp0_x = d->inter->interpol_d[arg[0]][i0];
	  for (i1=0; i1<p->cao; i1++) {
	    j = wrap_mesh_index(base[1] + i1, d->mesh);
	    tmp1 = d->inter->interpol[arg[1]][i1];
	    tmp1_y = d->inter->interpol_d[arg[1]][i1];
	    for (i2=0; i2<p->cao; i2++) {
	      k = wrap_mesh_index(base[2] + i2, d->mesh);
	      tmp2 = d->inter->interpol[arg[2]][i2];
	      tmp2_z = d->inter->interpol_d[arg[2]][i2];
	      cur_ca_frac_val = q * tmp0 * tmp1 * tmp2;
	      d->cf[ii][cf_cnt] = cur_ca_frac_val ;

	      d->dQdx[ii][cf_cnt] = tmp0_x * tmp1 * tmp2 * qLeni;
	      d->dQdy[ii][cf_cnt] = tmp0 * tmp1_y * tmp2 * qLeni;
	      d->dQdz[ii][cf_cnt] = tmp0 * tmp1 * tmp2_z * qLeni;

	      d->Qmesh[c_ind(i,j,k)+ii] += cur_ca_frac_val;
	      cf_cnt++;
	    }
	  }
        }
    }
}

void assign_charge_and_derivatives_real(system_t *s, parameters_t *p, data_t *d)
{
    int dim, i0, i1, i2;
    int id;
    FLOAT_TYPE tmp0, tmp1, tmp2;
    FLOAT_TYPE tmp0_x, tmp1_y, tmp2_z;
    /* position of a particle in local mesh units */
    FLOAT_TYPE pos;
    /* 1d-index of nearest mesh point */
    int nmp;
    /* index for caf interpolation grid */
    int arg[3];
    /* index, index jumps for rs_mesh array */
    FLOAT_TYPE cur_ca_frac_val;
    int cf_cnt;

    int base[3];
    int i,j,k;
    int Mesh = p->mesh;
    FLOAT_TYPE MI2 = 2.0*(FLOAT_TYPE)MaxInterpol;

    FLOAT_TYPE Hi = (double)Mesh/s->length;
    FLOAT_TYPE Leni = 1.0/s->length;
    FLOAT_TYPE q, qLeni;

    double pos_shift;

    pos_shift = (double)((p->cao-1)/2);
    /* particle position in mesh coordinates */
    for (id=0;id<s->nparticles;id++) {
        cf_cnt = id*p->cao3;
	q = s->q[id] ;
	qLeni = q * Leni;
        for (dim=0;dim<3;dim++) {
	  pos    = s->p->fields[dim][id]*Hi - pos_shift;
	  nmp = int_floor(pos + 0.5);
	  base[dim]  = wrap_mesh_index( nmp, d->mesh);
	  arg[dim] = int_floor((pos - nmp + 0.5)*MI2);
	  d->ca_ind[0][3*id + dim] = base[dim];
        }

        for (i0=0; i0<p->cao; i0++) {
	  i = wrap_mesh_index(base[0] + i0, d->mesh);
	  tmp0 = d->inter->interpol[arg[0]][i0];
	  tmp0_x = d->inter->interpol_d[arg[0]][i0];
	  for (i1=0; i1<p->cao; i1++) {
	    j = wrap_mesh_index(base[1] + i1, d->mesh);
	    tmp1 = d->inter->interpol[arg[1]][i1];
	    tmp1_y = d->inter->interpol_d[arg[1]][i1];
	    for (i2=0; i2<p->cao; i2++) {
	      k = wrap_mesh_index(base[2] + i2, d->mesh);
	      tmp2 = d->inter->interpol[arg[2]][i2];
	      tmp2_z = d->inter->interpol_d[arg[2]][i2];
	      cur_ca_frac_val = q * tmp0 * tmp1 * tmp2;
	      d->cf[0][cf_cnt] = cur_ca_frac_val ;

	      d->dQdx[0][cf_cnt] = tmp0_x * tmp1 * tmp2 * qLeni;
	      d->dQdy[0][cf_cnt] = tmp0 * tmp1_y * tmp2 * qLeni;
	      d->dQdz[0][cf_cnt] = tmp0 * tmp1 * tmp2_z * qLeni;

	      d->Qmesh[Mesh*(Mesh+2) * i + (Mesh+2) * j + k] += cur_ca_frac_val;
	      cf_cnt++;
	    }
	  }
        }
    }
}


// assign the forces obtained from k-space
void assign_forces_ad(double force_prefac, system_t *s, parameters_t *p, data_t *d, forces_t *f, int ii)
{
  int i,i0,i1,i2;
  int cf_cnt=0;
  int *base;
  int j,k,l;
  FLOAT_TYPE B;
  FLOAT_TYPE force_x, force_y, force_z;
  int cao = p->cao;

  cf_cnt=0;
  for (i=0; i<s->nparticles; i++) {
    force_x = force_y = force_z = 0.0;
    base = d->ca_ind[ii] + 3*i;
    cf_cnt = i*p->cao3;
    for (i0=0; i0<cao; i0++) {
      j = wrap_mesh_index(base[0] + i0, d->mesh);
      for (i1=0; i1<cao; i1++) {
	k = wrap_mesh_index(base[1] + i1, d->mesh);
	for (i2=0; i2<cao; i2++) {
	  l = wrap_mesh_index(base[2] + i2, d->mesh);

	  B = force_prefac*d->Qmesh[c_ind(j,k,l)+ii];
	  force_x -= B*d->dQdx[ii][cf_cnt];
	  force_y -= B*d->dQdy[ii][cf_cnt];
	  force_z -= B*d->dQdz[ii][cf_cnt];
	  cf_cnt++;
	}
      }
    }
    f->f_k->fields[0][i] += force_x;
    f->f_k->fields[1][i] += force_y;
    f->f_k->fields[2][i] += force_z;
    if (ii==1) {
      f->f_k->fields[0][i] *= 0.5;
      f->f_k->fields[1][i] *= 0.5;
      f->f_k->fields[2][i] *= 0.5;
    }
  }
}

// assign the forces obtained from k-space
void assign_forces_ad_real(double force_prefac, system_t *s, parameters_t *p, data_t *d, forces_t *f)
{
  int i,i0,i1,i2;
  int cf_cnt=0;
  int *base;
  int j,k,l;
  FLOAT_TYPE B;
  FLOAT_TYPE force_x, force_y, force_z;
  int cao = p->cao;
  int mesh = p->mesh;

  cf_cnt=0;
  for (i=0; i<s->nparticles; i++) {
    force_x = force_y = force_z = 0.0;
    base = d->ca_ind[0] + 3*i;
    cf_cnt = i*p->cao3;
    for (i0=0; i0<cao; i0++) {
      j = wrap_mesh_index(base[0] + i0, d->mesh);
      for (i1=0; i1<cao; i1++) {
	k = wrap_mesh_index(base[1] + i1, d->mesh);
	for (i2=0; i2<cao; i2++) {
	  l = wrap_mesh_index(base[2] + i2, d->mesh);

	  B = force_prefac*d->Qmesh[mesh*(mesh+2) * j + (mesh+2) * k + l];
	  force_x -= B*d->dQdx[0][cf_cnt];
	  force_y -= B*d->dQdy[0][cf_cnt];
	  force_z -= B*d->dQdz[0][cf_cnt];
	  cf_cnt++;
	}
      }
    }
    f->f_k->fields[0][i] += force_x;
    f->f_k->fields[1][i] += force_y;
    f->f_k->fields[2][i] += force_z;
  }
}

