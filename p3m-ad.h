#ifndef P3M_AD_H
#define P3M_AD_H

#include "p3m.h"

void Influence_function_berechnen_ad(FLOAT_TYPE alpha);
void P3M_ad(const FLOAT_TYPE alpha, const int Teilchenzahl);
void Init_ad(int Teilchenzahl);

const method_t method_p3m_ad = { 2, "P3M with analytic differentiation, not intelaced.", P3M_FLAG_ad, &Init_ad, &Influence_function_berechnen_ad, &P3M_ad, NULL };

#endif
