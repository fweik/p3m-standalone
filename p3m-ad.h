#pragma once

#ifndef P3M_AD_H
#define P3M_AD_H

#include "types.h"

void Influence_function_berechnen_ad( system_t *, parameters_t *, data_t * );
void P3M_ad( system_t *, parameters_t *, data_t *, forces_t * );
data_t *Init_ad( system_t *, parameters_t * );
FLOAT_TYPE Error_ad( system_t *, parameters_t * );

extern const method_t method_p3m_ad;

#endif
