#ifndef GENERATE_SYSTEM_H
#define GENERATE_SYSTEM_H

#include "common.h"

enum { FORM_FACTOR_RANDOM = 0,
       FORM_FACTOR_INNER_BOX = 1,
       FORM_FACTOR_MADELUNG = 2
     };

system_t *generate_system( int, int, FLOAT_TYPE, FLOAT_TYPE);

#endif
