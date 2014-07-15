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

#ifndef IO_H
#define IO_H

#include "types.h"

void Read_exact_forces(system_t *s, char *);
system_t *Read_system(parameters_t *, char *);
void Write_exact_forces(system_t *, char *);
void Write_system(system_t *, char *);
void Write_system_cuda( system_t *s, parameters_t *p, char *filename);

void write_vtf(char *filename, system_t *s);

void print_parameters(parameters_t p);

#endif
