#pragma once

#ifndef P3M_AD_I_H
#define P3M_AD_I_H

#include "types.h"

void Influence_function_ad_i( system_t *, parameters_t *, data_t * );
void P3M_ad_i( system_t *, parameters_t *, data_t *, forces_t * );
data_t *Init_ad_i( system_t *, parameters_t * );
FLOAT_TYPE Error_ad_i( system_t *, parameters_t * );

extern const method_t method_p3m_ad_i;

//const method_t method_p3m_ad_i = { METHOD_P3M_ad_i, "P3M with analytic differentiation, not intelaced.", METHOD_FLAG_ad | METHOD_FLAG_interlaced, &Init_ad_i, &Influence_function_berechnen_ad_i, &P3M_ad_i, NULL };

#endif
