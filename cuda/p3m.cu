#define PI 3.14159

#define P3M_BRILLOUIN 1
#define SQR(A) ((A)*(A))

__device__ double sinc(double d)
{
  double PId = PI*d;
  return (d == 0.0) ? 1.0 : sin(PId)/PId;
}

__global__ void influence_function( double *G_hat, double box, int cao, int mesh, double alpha ) {
  int n[3];
  int linear_index;
  double nom[3] = { 0.0, 0.0, 0.0 }, dnom = 0.0;
  double fak = SQR( PI / alpha );
  int mx, my, mz;
  double box_i = 1./box;
  int nshift[3];
  int nmx, nmy, nmz, nm2;
  double zwi;
  double S1, S2, S3;
  
  n[0] = blockDim.x * blockIdx.x;
  n[1] = blockDim.y * blockIdx.y;
  n[2] = threadIdx.x;

  nshift[0] = n[0] - round(n[0]/(double)mesh) * (double)mesh;
  nshift[1] = n[1] - round(n[1]/(double)mesh) * (double)mesh;
  nshift[2] = n[2] - round(n[2]/(double)mesh) * (double)mesh;

  linear_index = SQR(mesh)*n[0] + mesh * n[1] + n[2];

  if( (n[0] % (mesh/2))  && (n[1] % (mesh/2))  && (n[2] % (mesh/2))) {
    G_hat[linear_index] = 0.0;
    return;
  } 

  for ( mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++ ) {
    nmx = nshift[n[0]] + mesh*mx;
    S1 = pow ( sinc ( box_i*nmx ), 2*cao );
    for ( my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++ ) {
      nmy = nshift[n[1]] + mesh*my;
      S2   = S1*pow ( sinc ( box_i*nmy ), 2*cao );
      for ( mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++ ) {
	nmz = nshift[n[2]] + mesh*mz;
	S3   = S2*pow ( sinc ( box_i*nmz ), 2*cao );

	nm2 = SQR ( nmx*box_i ) + SQR ( nmy*box_i ) + SQR ( nmz*box_i );
	dnom += S3;

	zwi  = S3 * exp ( -fak*nm2 )/nm2;
	nom[0] += nmx*zwi*box_i;
	nom[1] += nmy*zwi*box_i;
	nom[2] += nmz*zwi*box_i;
      }
    }
  }
  
  zwi = box_i * (nshift[0]*nom[0] + nshift[1]*nom[1] + nshift[2]*nom[2]);
  zwi /= (SQR(nshift[0]) + SQR(nshift[1]) + SQR(nshift[2])) * SQR(box_i) *SQR(dnom);
  
  G_hat[linear_index] = 2.0 * zwi / PI;

  return;
}

__global__ void convolute( double *G_hat, double *mesh, int mesh_size ) {
  int n[3];
  int linear_index;
  n[0] = blockDim.x * blockIdx.x;
  n[1] = blockDim.y * blockIdx.y;
  n[2] = threadIdx.x;
  
  linear_index = SQR(mesh_size)*n[0] + mesh_size * n[1] + n[2];

  mesh[linear_index] *= G_hat[linear_index];
}

