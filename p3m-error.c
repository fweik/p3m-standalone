#include <math.h>

#include "p3m-error.h"

static double sinc(double d)
{
#define epsi 0.1

#define c2 -0.1666666666667e-0
#define c4  0.8333333333333e-2
#define c6 -0.1984126984127e-3
#define c8  0.2755731922399e-5

  double PId = PI*d, PId2;

  if (fabs(d)>epsi)
    return sin(PId)/PId;
  else {
    PId2 = SQR(PId);
    return 1.0 + PId2*(c2+PId2*(c4+PId2*(c6+PId2*c8)));
  }
  return 1.0;
}

static double analytic_cotangent_sum(int n, double mesh_i, int cao)
{
  double c, res=0.0;
  c = SQR(cos(PI*mesh_i*(double)n));

  switch (cao) {
  case 1 : { 
    res = 1; 
    break; }
  case 2 : { 
    res = (1.0+c*2.0)/3.0; 
    break; }
  case 3 : { 
    res = (2.0+c*(11.0+c*2.0))/15.0; 
    break; }
  case 4 : { 
    res = (17.0+c*(180.0+c*(114.0+c*4.0)))/315.0; 
    break; }
  case 5 : { 
    res = (62.0+c*(1072.0+c*(1452.0+c*(247.0+c*2.0))))/2835.0; 
    break; }
  case 6 : { 
    res = (1382.0+c*(35396.0+c*(83021.0+c*(34096.0+c*(2026.0+c*4.0)))))/155925.0; 
    break; }
  case 7 : { 
    res = (21844.0+c*(776661.0+c*(2801040.0+c*(2123860.0+c*(349500.0+c*(8166.0+c*4.0))))))/6081075.0; 
    break; }
  }
  
  return res;
}

void p3m_tune_aliasing_sums_ik(int nx, int ny, int nz, 
			    int mesh[3], double mesh_i[3], int cao, double alpha_L_i, 
			    double *alias1, double *alias2)
{

  int    mx,my,mz;
  double nmx,nmy,nmz;
  double fnmx,fnmy,fnmz;

  double ex,ex2,nm2,U2,factor1;

  factor1 = SQR(PI*alpha_L_i);

  *alias1 = *alias2 = 0.0;
  for (mx=-P3M_BRILLOUIN_TUNING; mx<=P3M_BRILLOUIN_TUNING; mx++) {
    fnmx = mesh_i[0] * (nmx = nx + mx*mesh[0]);
    for (my=-P3M_BRILLOUIN_TUNING; my<=P3M_BRILLOUIN_TUNING; my++) {
      fnmy = mesh_i[1] * (nmy = ny + my*mesh[1]);
      for (mz=-P3M_BRILLOUIN_TUNING; mz<=P3M_BRILLOUIN_TUNING; mz++) {
	fnmz = mesh_i[2] * (nmz = nz + mz*mesh[2]);

	nm2 = SQR(nmx) + SQR(nmy) + SQR(nmz);
	ex2 = SQR( ex = exp(-factor1*nm2) );
	
	U2 = pow(sinc(fnmx)*sinc(fnmy)*sinc(fnmz), 2.0*cao);
	
	*alias1 += ex2 / nm2;
	*alias2 += U2 * ex * (nx*nmx + ny*nmy + nz*nmz) / nm2;
      }
    }
  }
}

double p3m_real_space_error(double prefac, double r_cut_iL, 
			    int n_c_part, double sum_q2, double alpha_L,  double *box_l)
{
  return (2.0*prefac*sum_q2*exp(-SQR(r_cut_iL*alpha_L))) / (sqrt((double)n_c_part*r_cut_iL)*box_l[1]*box_l[2]);
}

double p3m_k_space_error_ik(double prefac, int mesh[3], int cao, int n_c_part, double sum_q2, double alpha_L, double *box_l)
{
  int  nx, ny, nz;
  double he_q = 0.0, mesh_i[3] = {1.0/mesh[0], 1.0/mesh[1], 1.0/mesh[2]}, alpha_L_i = 1./alpha_L;
  double alias1, alias2, n2, cs;
  double ctan_x, ctan_y;

  for (nx=-mesh[0]/2; nx<mesh[0]/2; nx++) {
    ctan_x = analytic_cotangent_sum(nx,mesh_i[0],cao);
    for (ny=-mesh[1]/2; ny<mesh[1]/2; ny++) {
      ctan_y = ctan_x * analytic_cotangent_sum(ny,mesh_i[1],cao);
      for (nz=-mesh[2]/2; nz<mesh[2]/2; nz++) {
	if((nx!=0) || (ny!=0) || (nz!=0)) {
	  n2 = SQR(nx) + SQR(ny) + SQR(nz);
	  cs = analytic_cotangent_sum(nz,mesh_i[2],cao)*ctan_y;
	  p3m_tune_aliasing_sums_ik(nx,ny,nz,mesh,mesh_i,cao,alpha_L_i,&alias1,&alias2);
	  he_q += (alias1  -  SQR(alias2/cs) / n2);
	}
      }
    }
  }
  return 2.0*prefac*sum_q2*sqrt(he_q/(double)n_c_part) / (box_l[1]*box_l[2]);
}

void p3m_tune_aliasing_sums_ad(int nx, int ny, int nz, 
			    int mesh, double mesh_i, int cao, double alpha_L_i, 
			    double *alias1, double *alias2, double *alias3,double *alias4)
{

  int    mx,my,mz;
  double nmx,nmy,nmz;
  double fnmx,fnmy,fnmz;

  double ex,ex2,nm2,U2,factor1;

  factor1 = SQR(PI*alpha_L_i);

  *alias1 = *alias2 = *alias3 = *alias4 = 0.0;
  for (mx=-P3M_BRILLOUIN_TUNING; mx<=P3M_BRILLOUIN_TUNING; mx++) {
    fnmx = mesh_i * (nmx = nx + mx*mesh);
    for (my=-P3M_BRILLOUIN_TUNING; my<=P3M_BRILLOUIN_TUNING; my++) {
      fnmy = mesh_i * (nmy = ny + my*mesh);
      for (mz=-P3M_BRILLOUIN_TUNING; mz<=P3M_BRILLOUIN_TUNING; mz++) {
	fnmz = mesh_i * (nmz = nz + mz*mesh);
	
	nm2 = SQR(nmx) + SQR(nmy) + SQR(nmz);
	ex2 = SQR( ex = exp(-factor1*nm2) );
	
	U2 = pow(sinc(fnmx)*sinc(fnmy)*sinc(fnmz), 2.0*cao);
	
	*alias1 += ex2 / nm2;
	*alias2 += U2 * ex;
	*alias3 += U2 * nm2;
	*alias4 += U2;
      }
    }
  }
}


double P3M_k_space_error_AD(double box_size, double prefac, int mesh, 
			 int cao, int n_c_part, double sum_q2, double alpha_L)
{
  int  nx, ny, nz;
  double he_q = 0.0, mesh_i = 1./mesh, alpha_L_i = 1./alpha_L;
  double alias1, alias2, alias3, alias4, n2, cs;

  for (nx=-mesh/2; nx<mesh/2; nx++)
    for (ny=-mesh/2; ny<mesh/2; ny++)
      for (nz=-mesh/2; nz<mesh/2; nz++)
	if((nx!=0) || (ny!=0) || (nz!=0)) {
	  n2 = SQR(nx) + SQR(ny) + SQR(nz);
	  p3m_tune_aliasing_sums_ad(nx,ny,nz,mesh,mesh_i,cao,alpha_L_i,&alias1,&alias2,&alias3,&alias4);	//alias4 = cs
	  he_q += (alias1  -  SQR(alias2) / (alias3*alias4));
	}
  return 2.0*prefac*sum_q2*sqrt(he_q/(double)n_c_part) / SQR(box_size);
}


double p3m_error_ik(double prefac, int mesh[3], int cao, int n_c_part, double sum_q2, double alpha_L, double r_cut_iL, double *box_l) {
  return sqrt(SQR(p3m_real_space_error(prefac, r_cut_iL,  n_c_part, sum_q2, alpha_L,  box_l)) \
    + SQR(p3m_k_space_error_ik(prefac, mesh, cao, n_c_part, sum_q2, alpha_L, box_l)));
}
