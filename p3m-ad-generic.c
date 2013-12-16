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

#include "types.h"
#include "p3m-common.h"

p3m_method_t *Init_p3m_method( R3_to_R A, R3_to_R B, int flags) {
  p3m_method_t *m = Init_array( 1, sizeof(p3m_method_t));

  m->A = A;
  m->B = B;

  return m;
}

void Init_influence_function_generic( p3m_metod_t *m, system_t *s, parameters_t *p, data_t *d) {
  int nx, ny, nz;

  for (nx=-d->mesh/2+1; nx<d->mesh/2; nx++) {
    for (ny=-d->mesh/2+1; ny<d->mesh/2; ny++) {
      for (nz=-d->mesh/2+1; nz<d->mesh/2; nz++) {
	if((nx!=0) && (ny!=0) && (nz!=0)) {
	  ind = r_ind(NTRANS(nx), NTRANS(ny), NTRANS(nz));
	  d->G_hat[ind] = m->B(nx,ny,nz,s,p) / m->A(nx,ny,nz,s,p);
	} else {
	  d->G_hat[ind] = 0.0;
	}
      }
    }
  }
}




data_t *Init_p3m_ad_generic( R3_to_R A, R3_to_R B, R3_to_R C, int flags) {

}
