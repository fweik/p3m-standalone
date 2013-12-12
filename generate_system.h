#ifndef GENERATE_SYSTEM_H
#define GENERATE_SYSTEM_H

#include "common.h"

enum { SYSTEM_RANDOM = 0,
       SYSTEM_INNER_BOX = 1,
       SYSTEM_MADELUNG = 2,
       SYSTEM_SEPARATED_DIPOLE = 3,
       SYSTEM_GAUSSIAN = 4,
       SYSTEM_SLAB = 5
     };

system_t *generate_system( int, int, FLOAT_TYPE, FLOAT_TYPE);

#endif
