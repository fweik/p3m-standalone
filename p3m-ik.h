#ifndef P3M_IK_H
#define P3M_IK_H

#define P3M_BRILLOUIN_TUNING 0

void Influence_function_berechnen_ik(FLOAT_TYPE alpha);
void P3M_ik(const FLOAT_TYPE alpha, const int Teilchenzahl);
void Init_ik(int Teilchenzahl);


const method_t method_p3m_ik = { 0, "P3M with ik differentiation, not intelaced.", P3M_FLAG_ik, &Init_ik, &Influence_function_berechnen_ik, &P3M_ik, NULL };

#endif
