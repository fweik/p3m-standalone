#pragma once

#ifndef P3M_IK_I_H
#define P3M_IK_I_H

#include "common.h"

void Influence_ik_i( system_t *, parameters_t *, data_t * );
void P3M_ik_i( system_t *, parameters_t *, data_t *, forces_t * );
data_t *Init_ik_i( system_t *, parameters_t * );
FLOAT_TYPE Error_ik_i( system_t *, parameters_t *);

extern const method_t method_p3m_ik_i;

#endif 
