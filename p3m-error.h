#ifndef P3M_ERROR_H
#define P3M_ERROR_H

#define P3M_BRILLOUIN_TUNING 0
#define PI 3.14159265358979323846264

#ifndef SQR
  #define SQR(A) ((A) * (A))
#endif

double p3m_error_ik(double, int*, int, int, double, double, double, double *);
double p3m_error_ad(double, int*, int, int, double, double, double, double *);

#endif
