#ifndef P3M_IK_H
#define P3M_IK_H

#include "p3m.h"

#define P3M_BRILLOUIN_TUNING 0

void Influence_function_berechnen_ik(FLOAT_TYPE alpha);
void P3M_ik(const FLOAT_TYPE alpha, const int Teilchenzahl);
void Init_ik(int Teilchenzahl);


extern const method_t method_p3m_ik;

#endif
