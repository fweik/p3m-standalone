#pragma once

#ifndef P3M_IK_I_H
#define P3M_IK_I_H

#include "common.h"

void Influence_function_berechnen_ik_i(FLOAT_TYPE alpha);
void P3M_ik_i(const FLOAT_TYPE alpha, const int Teilchenzahl);
void Init_ik_i(int Teilchenzahl);

const method_t method_p3m_ik_i = { METHOD_P3M_ik_i, "P3M with ik differentiation, intelaced.", METHOD_FLAG_ik | METHOD_FLAG_interlaced, &Init_ik_i, &Influence_function_berechnen_ik_i, &P3M_ik_i, NULL };

#endif 
