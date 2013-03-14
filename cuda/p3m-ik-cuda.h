#ifndef P3M_IK_CUDA_I_H
#define P3M_IK_CUDA_I_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int n;
  double *pos;
  double *q;
  double **f;
  double alpha;
  int cao;
  int mesh;
  double box;
  void *s;
} p3m_cuda_data_t;

void p3m_ik_cuda_init( p3m_cuda_data_t *d );
void p3m_ik_cuda( p3m_cuda_data_t *d );
void p3m_ik_cuda_free(void);
#ifdef __cplusplus
}
#endif
#endif
