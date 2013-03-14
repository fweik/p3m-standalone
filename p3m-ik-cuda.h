#ifndef P3M_IK_CUDA_H
#define P3M_IK_CUDA_H

#include "types.h"

extern const method_t method_p3m_ik_cuda;

data_t *Init_ik_cuda( system_t *s, parameters_t *p );
void P3M_ik_cuda( system_t *s, parameters_t *p, data_t *d, forces_t *f);


#endif
