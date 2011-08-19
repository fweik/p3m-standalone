#pragma once

#ifndef P3M_AD_I_H
#define P3M_AD_I_H

#include "p3m.h"

void Influence_function_berechnen_ad_i(FLOAT_TYPE);
void P3M_ad_i(FLOAT_TYPE,int);
void Init_ad_i(int);

const method_t method_p3m_ad_i = { METHOD_P3M_ad_i, "P3M with analytic differentiation, not intelaced.", P3M_FLAG_ad | P3M_FLAG_interlaced, &Init_ad_i, &Influence_function_berechnen_ad_i, &P3M_ad_i, NULL };

#endif
