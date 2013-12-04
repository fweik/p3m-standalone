#pragma once

#ifndef P3M_AD_R_H
#define P3M_AD_R_H

#include "types.h"

#define P3M_AD_SELF_FORCES

void Influence_function_berechnen_ad_r( system_t *, parameters_t *, data_t * );
void P3M_ad_r( system_t *, parameters_t *, data_t *, forces_t * );
data_t *Init_ad_r( system_t *, parameters_t * );
FLOAT_TYPE Error_ad_r( system_t *, parameters_t * );
FLOAT_TYPE p3m_k_space_error_ad_r( system_t *, parameters_t * );

extern const method_t method_p3m_ad_r;

#endif
