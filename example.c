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

// Methods

#include "p3m-ik.h"

// Utils and IO

#include "io.h"

// Error calculation

#include "error.h"

// Real space methods

#include "realpart.h"

int main ( void ) {
    system_t *system;
    parameters_t parameters;
    data_t *data;
    forces_t *forces;
    error_t error;

    // Inits the system and reads particle data and parameters from file.
    system = Read_system ( &parameters, "tests-wall.pos" );
    forces = Init_forces(system->nparticles);

    // Calculate the reference forces
    Calculate_reference_forces( system, &parameters );

    // Init data structures for the method
    data = method_p3m_ik.Init ( system, &parameters );

    // Init the neighbor list for the real space part
    Init_neighborlist ( system, &parameters, data );

    // Calculate the forces
    Calculate_forces ( &method_p3m_ik, system, &parameters, data, forces ); /* Hockney/Eastwood */

    // Calculate the error
    error = Calculate_errors ( system, forces );

    return 0;
}

