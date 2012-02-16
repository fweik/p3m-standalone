#include "window-functions.h"

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

FLOAT_TYPE caf_bspline_k(int i, FLOAT_TYPE d)
{
#define epsi 0.1

#define c2 -0.1666666666667e-0
#define c4  0.8333333333333e-2
#define c6 -0.1984126984127e-3
#define c8  0.2755731922399e-5

    double PId = PI*d, PId2;

    if (fabs(d)>epsi)
      return pow(sin(PId)/PId, i);
    else {
        PId2 = SQR(PId);
        return pow(1.0 + PId2*(c2+PId2*(c4+PId2*(c6+PId2*c8))), i);
    }
    return 1.0;
}

FLOAT_TYPE caf_bspline_d(int i, FLOAT_TYPE x, int cao_value) {
  int ip = cao_value - 1;
  switch (ip) 
    {
    case 1:
	    switch(i) {
	    case 0: return -1.0;
	    case 1: return 1.0;
	    }
      break;
    case 2:
	    switch(i) {
	    case 0: return x-0.5;
	    case 1: return -2.0*x;
	    case 2: return x+0.5;
	    }
      break;
    case 3:
	    switch(i) {
	    case 0: return (-1.0+x*(  4.0+x*( -4.0)))/8.0;
	    case 1: return (-5.0+x*( -4.0+x*( 12.0)))/8.0;
	    case 2: return ( 5.0+x*( -4.0+x*(-12.0)))/8.0;
	    case 3: return ( 1.0+x*(  4.0+x*(  4.0)))/8.0;
	    }
      break;
    case 4:
	    switch(i) {
	    case 0: return ( -1.0+x*(  6.0+x*( -12.0+x*( 8.0))))/48.0;
	    case 1: return (-11.0+x*( 12.0+x*(  12.0+x*(-16.0))))/24.0;
	    case 2: return (      x*(-5.0+x*x*4.0))/4.0;
	    case 3: return ( 11.0+x*( 12.0+x*( -12.0+x*(-16.0))))/ 24.0;
	    case 4: return (  1.0+x*(  6.0+x*(  12.0+x*(  8.0))))/48.0;
	    }
        break;
    case 5:
	    switch(i) {
	    case 0: return ( -1.0+x*(  8.0+x*( -24.0+x*(  32.0+x*(-16)))))/384.0;
	    case 1: return (-75.0+x*( 168.0+x*( -72.0+x*( -96.0+x*(80.0)))))/384.0;
	    case 2: return (-77.0+x*( -88.0+x*( 168.0+x*(  32.0+x*(-80.0)))))/192.0;
	    case 3: return ( 77.0+x*( -88.0+x*(-168.0+x*(  32.0+x*(80.0)))))/192.0;
	    case 4: return ( 75.0+x*( 168.0+x*(  72.0+x*( -96.0+x*(-80)))))/384.0;
	    case 5: return (  1.0+x*(   8.0+x*(  24.0+x*(  32.0+x*(16.0)))))/384.0;
	    }
       break;
    case 6 : 
	    switch(i) {
	    case 0: return (  -1.0+x*( 10.0+x*( -40.0+x*(  80.0+x*(-80.0+x*32.0)))))/3840.0;
	    case 1: return ( -59.0+x*(185.0+x*(-200.0+x*(  40.0+x*( 80.0-x*48.0)))))/960.0;
	    case 2: return (-289.0+x*(158.0+x*( 344.0+x*(-272.0+x*(-80.0+x*96.0)))))/768.0;
	    case 3: return (       x*(-77.0+        x*x*(  56.0       -x*x*16.0) ) )/96.0;
	    case 4: return ( 289.0+x*(158.0+x*(-344.0+x*(-272.0+x*( 80.0+x*96.0)))))/768.0;
	    case 5: return (  59.0+x*(185.0+x*( 200.0+x*(  40.0+x*(-80.0-x*48.0)))))/960.0;
	    case 6: return (   1.0+x*( 10.0+x*(  40.0+x*(  80.0+x*( 80.0+x*32.0)))))/3840.0;
	    }
       break;
    default :
      fprintf(stderr,"Charge assignment order %d unknown.\n",cao_value);
      return 0.0;
  }
  return 0.0;
}

FLOAT_TYPE caf_bspline(int i, FLOAT_TYPE x, int cao_value) {
  switch (cao_value) {
  case 1 : return 1.0;
  case 2 : {
    switch (i) {
    case 0: return 0.5-x;
    case 1: return 0.5+x;
    default:
      fprintf(stderr,"Tried to access charge assignment function of degree %d in scheme of order %d.\n",i,cao_value);
      return 0.0;
    }
  } 
  case 3 : { 
    switch (i) {
    case 0: return 0.5*SQR(0.5 - x);
    case 1: return 0.75 - SQR(x);
    case 2: return 0.5*SQR(0.5 + x);
    default:
      fprintf(stderr,"Tried to access charge assignment function of degree %d in scheme of order %d.\n",i,cao_value);
      return 0.0;
    }
  case 4 : { 
    switch (i) {
    case 0: return ( 1.0+x*( -6.0+x*( 12.0-x* 8.0)))/48.0;
    case 1: return (23.0+x*(-30.0+x*(-12.0+x*24.0)))/48.0;
    case 2: return (23.0+x*( 30.0+x*(-12.0-x*24.0)))/48.0;
    case 3: return ( 1.0+x*(  6.0+x*( 12.0+x* 8.0)))/48.0;
    default:
      fprintf(stderr,"Tried to access charge assignment function of degree %d in scheme of order %d.\n",i,cao_value);
      return 0.0;
    }
  }
  case 5 : {
    switch (i) {
    case 0: return (  1.0+x*( -8.0+x*(  24.0+x*(-32.0+x*16.0))))/384.0;
    case 1: return ( 19.0+x*(-44.0+x*(  24.0+x*( 16.0-x*16.0))))/ 96.0;
    case 2: return (115.0+x*       x*(-120.0+x*       x*48.0))  /192.0;
    case 3: return ( 19.0+x*( 44.0+x*(  24.0+x*(-16.0-x*16.0))))/ 96.0;
    case 4: return (  1.0+x*(  8.0+x*(  24.0+x*( 32.0+x*16.0))))/384.0;
    default:
      fprintf(stderr,"Tried to access charge assignment function of degree %d in scheme of order %d.\n",i,cao_value);
      return 0.0;
    }
  }
  case 6 : {
    switch (i) {
    case 0: return (  1.0+x*( -10.0+x*(  40.0+x*( -80.0+x*(  80.0-x* 32.0)))))/3840.0;
    case 1: return (237.0+x*(-750.0+x*( 840.0+x*(-240.0+x*(-240.0+x*160.0)))))/3840.0;
    case 2: return (841.0+x*(-770.0+x*(-440.0+x*( 560.0+x*(  80.0-x*160.0)))))/1920.0;
    case 3: return (841.0+x*(+770.0+x*(-440.0+x*(-560.0+x*(  80.0+x*160.0)))))/1920.0;
    case 4: return (237.0+x*( 750.0+x*( 840.0+x*( 240.0+x*(-240.0-x*160.0)))))/3840.0;
    case 5: return (  1.0+x*(  10.0+x*(  40.0+x*(  80.0+x*(  80.0+x* 32.0)))))/3840.0;
    default:
      fprintf(stderr,"Tried to access charge assignment function of degree %d in scheme of order %d.\n",i,cao_value);
      return 0.0;
    }
  }
  case 7 : {
    switch (i) {
    case 0: return (    1.0+x*(   -12.0+x*(   60.0+x*( -160.0+x*(  240.0+x*(-192.0+x* 64.0))))))/46080.0;
    case 1: return (  361.0+x*( -1416.0+x*( 2220.0+x*(-1600.0+x*(  240.0+x*( 384.0-x*192.0))))))/23040.0;
    case 2: return (10543.0+x*(-17340.0+x*( 4740.0+x*( 6880.0+x*(-4080.0+x*(-960.0+x*960.0))))))/46080.0;
    case 3: return ( 5887.0+x*          x*(-4620.0+x*         x*( 1680.0-x*        x*320.0)))   /11520.0;
    case 4: return (10543.0+x*( 17340.0+x*( 4740.0+x*(-6880.0+x*(-4080.0+x*( 960.0+x*960.0))))))/46080.0;
    case 5: return (  361.0+x*(  1416.0+x*( 2220.0+x*( 1600.0+x*(  240.0+x*(-384.0-x*192.0))))))/23040.0;
    case 6: return (    1.0+x*(    12.0+x*(   60.0+x*(  160.0+x*(  240.0+x*( 192.0+x* 64.0))))))/46080.0;
    default:
      fprintf(stderr,"Tried to access charge assignment function of degree %d in scheme of order %d.\n",i,cao_value);
      return 0.0;
    }
  }
  default :{
    fprintf(stderr,"Charge assignment order %d unknown.\n",cao_value);
    return 0.0;
  }}}
}
