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

#ifndef P3M_AD_SELF_FORCES_H
#define P3M_AD_SELF_FORCES_H

#include "types.h"

//#define P3M_AD_SELF_FORCES_DEBUG

#ifdef P3M_AD_SELF_FORCES_DEBUG
#define SF_TRACE(A); A;
#else
#define SF_TRACE(A)
#endif


void Init_self_forces( system_t *s, parameters_t *p, data_t *d );
void Substract_self_forces( system_t *s, parameters_t *p, data_t *d, forces_t *f );

#endif
