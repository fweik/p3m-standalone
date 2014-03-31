#include <math.h>
#include <stdio.h>

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062
#define SQR(A) ((A)*(A))

double sinc(double x) {
  return (x==0)?1.0:sin(PI*x)/(PI*x);
}

double U(double x, double y, double z int p) {
  return pow(sinc(x)*sinc(y)*sinc(z), p);
}

FLOAT_TYPE K2(int nx, int ny, int nz, FLOAT_TYPE l) {
  return SQR(2.0*PI/l) * ( SQR ( nx ) + SQR ( ny ) + SQR ( nz ) ); 
}

FLOAT_TYPE Phi(int nx, int ny, int nz, FLOAT_TYPE l, FLOAT_TYPE alpha) {
  FLOAT_TYPE k2 = K2(nx, ny, nz, l);
  return (4*PI/k2) * exp(-k2/(4*alpha*alpha));
}

FLOAT_TYPE C(int nx, int ny, int nz, FLOAT_TYPE l, FLOAT_TYPE alpha) {
  FLOAT_TYPE k2 = K2(nx, ny, nz,l);
  FLOAT_TYPE phi = Phi(nx,ny,nz,l,alpha);

  return k2 * SQR(phi);
}


int main() {
  int i,j,k;
  int mesh=32, cao=5;
  double l=10.0;
  double alpha;
  int mc=1;
  int mx,my,mz;
  double s1=0.0, s2 = 0.0;
  double G_trunc;
  double a1, a2, u2;

  for(alpha=0.0;alpha<=3.0;alpha+=0.01) {
    for(i=-mesh/2;i<=mesh/2;i++)
      for(j=-mesh/2;j<=mesh/2;j++)
	for(k=-mesh/2;k<=mesh/2;k++)  {
	  a1 = a2 = 0.0;
	  for(mx=-mc;mx<=mc;mx++) 
	    for(mx=-mc;mx<=mc;mx++) 
	      for(mx=-mc;mx<=mc;mx++) {
		u2 = U(i/mesh+mc,j/mesh+my, k/mesh+mz);
		a1 += u2;
		a2 += K2(i+mesh*mx, j+mesh*my, k+mesh*mz) * u2;
	      }
	}
  }
}
