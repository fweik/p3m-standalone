#pragma once

#ifndef P3M_AD_H
#define P3M_AD_H

#include "types.h"


void Influence_function_berechnen_ad( system_t *, parameters_t *, data_t * );
void P3M_ad( system_t *, parameters_t *, data_t *, forces_t * );
data_t *Init_ad( system_t *, parameters_t * );
FLOAT_TYPE Error_ad( system_t *, parameters_t * );
FLOAT_TYPE p3m_k_space_error_ad( system_t *, parameters_t * );

// Coefficients for the error-function
FLOAT_TYPE A_ad(int nx, int ny, int nz, system_t *s, parameters_t *p);
FLOAT_TYPE B_ad(int nx, int ny, int nz, system_t *s, parameters_t *p);

FLOAT_TYPE A_ad_dip(int nx, int ny, int nz, system_t *s, parameters_t *p);
FLOAT_TYPE B_ad_dip(int nx, int ny, int nz, system_t *s, parameters_t *p);

extern const method_t method_p3m_ad;

#endif
