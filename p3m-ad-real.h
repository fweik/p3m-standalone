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

#pragma once

#ifndef P3M_AD_R_H
#define P3M_AD_R_H

#include "types.h"

//#define P3M_AD_SELF_FORCES

void Influence_function_berechnen_ad_r( system_t *, parameters_t *, data_t * );
void P3M_ad_r( system_t *, parameters_t *, data_t *, forces_t * );
data_t *Init_ad_r( system_t *, parameters_t * );
FLOAT_TYPE Error_ad_r( system_t *, parameters_t * );
FLOAT_TYPE p3m_k_space_error_ad_r( system_t *, parameters_t * );

extern const method_t method_p3m_ad_r;

#endif
