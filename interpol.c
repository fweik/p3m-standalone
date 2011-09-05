#include <stdio.h>
#include <stdlib.h>

#include "types.h"
#include "interpol.h"
#include "common.h"

void Init_interpolation(int ip, data_t *d)
{
  /* 
     Die Ordnung ist an Petersens ip angelehnt.
     Im H/E-Fall ist die fouriertransformierte
     charge assignment function gerade sinc^(ip+1).
     (s.a. H/E, S. 159 f.)
  */
  
  FLOAT_TYPE dInterpol=(FLOAT_TYPE)MaxInterpol;
  FLOAT_TYPE x;
  long   i;

  FLOAT_TYPE **LadInt = d->LadInt = Init_array( ip + 2, sizeof(FLOAT_TYPE *));
  FLOAT_TYPE **LadInt_ = d->LadInt_ = Init_array( ip + 2, sizeof(FLOAT_TYPE *));

  for( i = 0; i <= ip; i++) {
    LadInt[i] = Init_array( 2* MaxInterpol + 1, sizeof(FLOAT_TYPE) );
    LadInt_[i] = Init_array( 2* MaxInterpol + 1, sizeof(FLOAT_TYPE) );
  }
  LadInt[ip+1] = LadInt_[ip+1] = NULL;

  /* charge assignment function: */
  switch (ip) 
    {
    case 0 : 
      { 
	for (i=-MaxInterpol; i<=MaxInterpol; i++)
	  {
	    x=i/(2.0*dInterpol);
	    LadInt[0][i+MaxInterpol] = 1.0;
	  }
      } break;
    case 1 : 
      { 
	for (i=-MaxInterpol; i<=MaxInterpol; i++)
	  {
	    x=i/(2.0*dInterpol);
	    LadInt[0][i+MaxInterpol] = 0.5-x;
	    LadInt[1][i+MaxInterpol] = 0.5+x;
	  }
      } break;
    case 2 : 
      { 
	for (i=-MaxInterpol; i<=MaxInterpol; i++)
	  {
	    x=i/(2.0*dInterpol);
	    LadInt[0][i+MaxInterpol] = 0.5*SQR(0.5 - x);
	    LadInt[1][i+MaxInterpol] = 0.75 - SQR(x);
	    LadInt[2][i+MaxInterpol] = 0.5*SQR(0.5 + x);
	  }
      } break;
    case 3 :
      { 
	for (i=-MaxInterpol; i<=MaxInterpol; i++)
	  {
	    x=i/(2.0*dInterpol);
	    LadInt[0][i+MaxInterpol] = ( 1.0+x*( -6.0+x*( 12.0-x* 8.0)))/48.0;
	    LadInt[1][i+MaxInterpol] = (23.0+x*(-30.0+x*(-12.0+x*24.0)))/48.0;
	    LadInt[2][i+MaxInterpol] = (23.0+x*( 30.0+x*(-12.0-x*24.0)))/48.0;
	    LadInt[3][i+MaxInterpol] = ( 1.0+x*(  6.0+x*( 12.0+x* 8.0)))/48.0;
	  }
      } break;
    case 4 : 
      {
	for (i=-MaxInterpol; i<=MaxInterpol; i++)
	  {
	    x=i/(2.0*dInterpol);
	    LadInt[0][i+MaxInterpol] = (  1.0+x*( -8.0+x*(  24.0+x*(-32.0+x*16.0))))/384.0;
	    LadInt[1][i+MaxInterpol] = ( 19.0+x*(-44.0+x*(  24.0+x*( 16.0-x*16.0))))/ 96.0;
	    LadInt[2][i+MaxInterpol] = (115.0+x*       x*(-120.0+x*       x*48.0))  /192.0;
	    LadInt[3][i+MaxInterpol] = ( 19.0+x*( 44.0+x*(  24.0+x*(-16.0-x*16.0))))/ 96.0;
	    LadInt[4][i+MaxInterpol] = (  1.0+x*(  8.0+x*(  24.0+x*( 32.0+x*16.0))))/384.0;
	  }
      } break;
    case 5 : 
      {
	for (i=-MaxInterpol; i<=MaxInterpol; i++)
	  {
	    x=i/(2.0*dInterpol);
	    LadInt[0][i+MaxInterpol] = (  1.0+x*( -10.0+x*(  40.0+x*( -80.0+x*(  80.0-x* 32.0)))))/3840.0;
	    LadInt[1][i+MaxInterpol] = (237.0+x*(-750.0+x*( 840.0+x*(-240.0+x*(-240.0+x*160.0)))))/3840.0;
	    LadInt[2][i+MaxInterpol] = (841.0+x*(-770.0+x*(-440.0+x*( 560.0+x*(  80.0-x*160.0)))))/1920.0;
	    LadInt[3][i+MaxInterpol] = (841.0+x*(+770.0+x*(-440.0+x*(-560.0+x*(  80.0+x*160.0)))))/1920.0;
	    LadInt[4][i+MaxInterpol] = (237.0+x*( 750.0+x*( 840.0+x*( 240.0+x*(-240.0-x*160.0)))))/3840.0;
	    LadInt[5][i+MaxInterpol] = (  1.0+x*(  10.0+x*(  40.0+x*(  80.0+x*(  80.0+x* 32.0)))))/3840.0;
	  }
      } break;
    case 6 :
      {
	for (i=-MaxInterpol; i<=MaxInterpol; i++)
	  {
	    x=i/(2.0*dInterpol);
	    LadInt[0][i+MaxInterpol] = (    1.0+x*(   -12.0+x*(   60.0+x*( -160.0+x*(  240.0+x*(-192.0+x* 64.0))))))/46080.0;
	    LadInt[1][i+MaxInterpol] = (  361.0+x*( -1416.0+x*( 2220.0+x*(-1600.0+x*(  240.0+x*( 384.0-x*192.0))))))/23040.0;
	    LadInt[2][i+MaxInterpol] = (10543.0+x*(-17340.0+x*( 4740.0+x*( 6880.0+x*(-4080.0+x*(-960.0+x*960.0))))))/46080.0;
	    LadInt[3][i+MaxInterpol] = ( 5887.0+x*          x*(-4620.0+x*         x*( 1680.0-x*        x*320.0)))   /11520.0;
	    LadInt[4][i+MaxInterpol] = (10543.0+x*( 17340.0+x*( 4740.0+x*(-6880.0+x*(-4080.0+x*( 960.0+x*960.0))))))/46080.0;
	    LadInt[5][i+MaxInterpol] = (  361.0+x*(  1416.0+x*( 2220.0+x*( 1600.0+x*(  240.0+x*(-384.0-x*192.0))))))/23040.0;
	    LadInt[6][i+MaxInterpol] = (    1.0+x*(    12.0+x*(   60.0+x*(  160.0+x*(  240.0+x*( 192.0+x* 64.0))))))/46080.0;
	  }
      } break;
    default :
      {
	fprintf(stderr,"Das Programm kennt fuer das H/E charge-assignment-scheme\n");
	fprintf(stderr,"den ip-Wert %d nicht! Programm abgebrochen.\n\n",ip);
	exit(1);
      }
    }

  switch (ip) 
    {
    case 1:
 	for (i=-MaxInterpol; i<=MaxInterpol; i++)
	  {
	    x=i/(2.0*dInterpol);
	    LadInt_[0][i+MaxInterpol] = -1.0;
	    LadInt_[1][i+MaxInterpol] = 1.0;
	  }
      break;
    case 2:
    	for (i=-MaxInterpol; i<=MaxInterpol; i++)
	  {
	    x=i/(2.0*dInterpol);
	    LadInt_[0][i+MaxInterpol] = x-0.5;
	    LadInt_[1][i+MaxInterpol] = -2.0*x;
	    LadInt_[2][i+MaxInterpol] = x+0.5;
	  }
      break;
    case 3:
	for (i=-MaxInterpol; i<=MaxInterpol; i++)
	  {
	    x=i/(2.0*dInterpol);
	    LadInt_[0][i+MaxInterpol] = (-1.0+x*(  4.0+x*( -4.0)))/8.0;
	    LadInt_[1][i+MaxInterpol] = (-5.0+x*( -4.0+x*( 12.0)))/8.0;
	    LadInt_[2][i+MaxInterpol] = ( 5.0+x*( -4.0+x*(-12.0)))/8.0;
	    LadInt_[3][i+MaxInterpol] = ( 1.0+x*(  4.0+x*(  4.0)))/8.0;
	  }
      break;
    case 4:
    	for (i=-MaxInterpol; i<=MaxInterpol; i++)
	  {
	    x=i/(2.0*dInterpol);
	    LadInt_[0][i+MaxInterpol] = ( -1.0+x*(  6.0+x*( -12.0+x*( 8.0))))/48.0;
	    LadInt_[1][i+MaxInterpol] = (-11.0+x*( 12.0+x*(  12.0+x*(-16.0))))/24.0;
	    LadInt_[2][i+MaxInterpol] = (      x*(-5.0+x*x*4.0))/4.0;
	    LadInt_[3][i+MaxInterpol] = ( 11.0+x*( 12.0+x*( -12.0+x*(-16.0))))/ 24.0;
	    LadInt_[4][i+MaxInterpol] = (  1.0+x*(  6.0+x*(  12.0+x*(  8.0))))/48.0;
	  }
        break;
    case 5:
	for (i=-MaxInterpol; i<=MaxInterpol; i++)
	  {
	    x=i/(2.0*dInterpol);
	    LadInt_[0][i+MaxInterpol] = ( -1.0+x*(  8.0+x*( -24.0+x*(  32.0+x*(-16)))))/384.0;
	    LadInt_[1][i+MaxInterpol] = (-75.0+x*( 168.0+x*( -72.0+x*( -96.0+x*(80.0)))))/384.0;
	    LadInt_[2][i+MaxInterpol] = (-77.0+x*( -88.0+x*( 168.0+x*(  32.0+x*(-80.0)))))/192.0;
	    LadInt_[3][i+MaxInterpol] = ( 77.0+x*( -88.0+x*(-168.0+x*(  32.0+x*(80.0)))))/192.0;
	    LadInt_[4][i+MaxInterpol] = ( 75.0+x*( 168.0+x*(  72.0+x*( -96.0+x*(-80)))))/384.0;
	    LadInt_[5][i+MaxInterpol] = (  1.0+x*(   8.0+x*(  24.0+x*(  32.0+x*(16.0)))))/384.0;
	    /*
	    x=i/(2.0*dInterpol);
	    LadInt_[0][i+MaxInterpol] = ( -1.0+x*( 8.0+x*( -24.0+x*(  32.0+x*(-16.0)))))/384.0;
	    LadInt_[1][i+MaxInterpol] = (-75.0+x*(168.0+x*(-72.0+x*(-96.0+x*( 80.0)))))/384.0;
	    LadInt_[2][i+MaxInterpol] = (-77.0+x*(-88.0+x*( 168.0+x*(32.0+x*(-80.0)))))/192.0;
	    LadInt_[3][i+MaxInterpol] = (-77.0+x*(-88.0+x*(-168.0+x*(32.0+x*(80.0) ))))/192.0;
	    LadInt_[4][i+MaxInterpol] = ( 75.0+x*(168.0+x*(72.0+x*(-96.0+x*(-80.0)))))/384.0;
	    LadInt_[5][i+MaxInterpol] = (  1.0+x*(  8.0+x*(  24.0+x*( 32.0+x*(16.0)))))/384.0;
	    */
	  }
       break;
    case 6 : 
      { 
	for (i=-MaxInterpol; i<=MaxInterpol; i++)
	  {
	    x=i/(2.0*dInterpol);
	    LadInt_[0][i+MaxInterpol] = (  -1.0+x*( 10.0+x*( -40.0+x*(  80.0+x*(-80.0+x*32.0)))))/3840.0;
	    LadInt_[1][i+MaxInterpol] = ( -59.0+x*(185.0+x*(-200.0+x*(  40.0+x*( 80.0-x*48.0)))))/960.0;
	    LadInt_[2][i+MaxInterpol] = (-289.0+x*(158.0+x*( 344.0+x*(-272.0+x*(-80.0+x*96.0)))))/768.0;
	    LadInt_[3][i+MaxInterpol] = (       x*(-77.0+        x*x*(  56.0       -x*x*16.0) ) )/96.0;
	    LadInt_[4][i+MaxInterpol] = ( 289.0+x*(158.0+x*(-344.0+x*(-272.0+x*( 80.0+x*96.0)))))/768.0;
	    LadInt_[5][i+MaxInterpol] = (  59.0+x*(185.0+x*( 200.0+x*(  40.0+x*(-80.0-x*48.0)))))/960.0;
	    LadInt_[6][i+MaxInterpol] = (   1.0+x*( 10.0+x*(  40.0+x*(  80.0+x*( 80.0+x*32.0)))))/3840.0;
	  }
      } break;
    default :
      {
	fprintf(stderr,"Das Programm kennt fuer die ABLEITUNG des H/E charge\n");
	fprintf(stderr,"assignment schemes den ip-Wert %d nicht! Programm abgebrochen.\n\n",ip);
	exit(1);
      }
    }
}
