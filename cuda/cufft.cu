#include <iostream>
#include <fstream>

#include <cufft.h>
#include <stdio.h>

#define N 8

#define PI 3.14159265359

#define SQR(A) ((A)*(A))
#define P3M_BRILLOUIN 1

using namespace std;

typedef struct {
  int n;
  double *pos;
  double *q;
  double *f_x;
  double *f_y;
  double *f_z;
  double alpha;
  int cao;
  int mesh;
  double box;
} data_t;

data_t *read_reference( char *filename ) {
  ifstream f;
  int i=0;
  data_t *d = (data_t *)malloc(sizeof(data_t));

  f.open(filename);

  f >> d->n;
  f >> d->cao;
  f >> d->mesh;
  f >> d->alpha;
  f >> d->box;

  d->pos = (double *)malloc(3*d->n*sizeof(double));
  d->q = (double *)malloc(d->n*sizeof(double));
  d->f_x = (double *)malloc(d->n*sizeof(double));
  d->f_y = (double *)malloc(d->n*sizeof(double));
  d->f_z = (double *)malloc(d->n*sizeof(double));

  while(f.good()) {
    f >> d->pos[3*i + 0];
    f >> d->pos[3*i + 1];
    f >> d->pos[3*i + 2];
    f >> d->q[i];
    f >> d->f_x[i];
    f >> d->f_y[i];
    f >> d->f_z[i];
    i++;
  }
  if(i != d->n)
    printf("Warning, not enought particles in file. (%d of %d)\n", i, d->n);

  return d;
}


__device__ __host__ inline double sinc(double d)
{
  double PId = PI*d;
  return (d == 0.0) ? 1.0 : sin(PId)/PId;
}

void Aliasing_sums_ik ( int cao, double box, double alpha, int mesh, int NX, int NY, int NZ,
                        double *Zaehler, double *Nenner ) {
    double S1,S2,S3;
    double fak1,fak2,zwi;
    int    MX,MY,MZ;
    double NMX,NMY,NMZ;
    double NM2;
    double expo, TE;
    double Leni = 1.0/box;

    fak1 = 1.0/ ( double ) mesh;
    fak2 = SQR ( PI/ ( alpha ) );

    Zaehler[0] = Zaehler[1] = Zaehler[2] = *Nenner = 0.0;

    for ( MX = -P3M_BRILLOUIN; MX <= P3M_BRILLOUIN; MX++ ) {
      NMX = ( ( NX > mesh/2 ) ? NX - mesh : NX ) + mesh*MX;
      S1 = pow ( sinc(fak1*NMX ), 2*cao );
      for ( MY = -P3M_BRILLOUIN; MY <= P3M_BRILLOUIN; MY++ ) {
	NMY = ( ( NY > mesh/2 ) ? NY - mesh : NY ) + mesh*MY;
	S2   = S1*pow ( sinc (fak1*NMY ), 2*cao );
	for ( MZ = -P3M_BRILLOUIN; MZ <= P3M_BRILLOUIN; MZ++ ) {
	  NMZ = ( ( NZ > mesh/2 ) ? NZ - mesh : NZ ) + mesh*MZ;
	  S3   = S2*pow ( sinc( fak1*NMZ ), 2*cao );

	  NM2 = SQR ( NMX*Leni ) + SQR ( NMY*Leni ) + SQR ( NMZ*Leni );
	  *Nenner += S3;

	  expo = fak2*NM2;
	  TE = exp ( -expo );
	  zwi  = S3 * TE/NM2;
	  Zaehler[0] += NMX*zwi*Leni;
	  Zaehler[1] += NMY*zwi*Leni;
	  Zaehler[2] += NMZ*zwi*Leni;
	}
      }
    }
}

/* Calculate influence function */
void Influence_function_berechnen_ik ( int cao, int mesh, double box, double alpha, double *G_hat ) {

  int    NX,NY,NZ;
  double Dnx,Dny,Dnz;
  double Zaehler[3]={0.0,0.0,0.0},Nenner=0.0;
  double zwi;
  int ind = 0;
  double Leni = 1.0/box;

  for ( NX=0; NX<mesh; NX++ ) {
    for ( NY=0; NY<mesh; NY++ ) {
      for ( NZ=0; NZ<mesh; NZ++ ) {
	ind = NX*mesh*mesh + NY * mesh + NZ;
	  
	if ( ( NX==0 ) && ( NY==0 ) && ( NZ==0 ) )
	  G_hat[ind]=0.0;
	else if ( ( NX% ( mesh/2 ) == 0 ) && ( NY% ( mesh/2 ) == 0 ) && ( NZ% ( mesh/2 ) == 0 ) )
	  G_hat[ind]=0.0;
	else {
	  Aliasing_sums_ik ( cao, box, alpha, mesh, NX, NY, NZ, Zaehler, &Nenner );
		  
	  Dnx = ( NX > mesh/2 ) ? NX - mesh : NX;
	  Dny = ( NY > mesh/2 ) ? NY - mesh : NY;
	  Dnz = ( NZ > mesh/2 ) ? NZ - mesh : NZ;
	    
	  zwi  = Dnx*Zaehler[0]*Leni + Dny*Zaehler[1]*Leni + Dnz*Zaehler[2]*Leni;
	  zwi /= ( ( SQR ( Dnx*Leni ) + SQR ( Dny*Leni ) + SQR ( Dnz*Leni ) ) * SQR ( Nenner ) );
	  G_hat[ind] = 2.0 * zwi / PI;
	}
      }
    }
  }
}


__device__ inline int wrap_index(const int ind, const int mesh) {
  if(ind < 0)
    return ind + mesh;
  else if(ind >= mesh)
    return ind - mesh;
  else 
    return ind;	   
}

__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

__device__ float caf(int i, float x, int cao_value) {
  switch (cao_value) {
  case 1 : return 1.0;
  case 2 : {
    switch (i) {
    case 0: return 0.5-x;
    case 1: return 0.5+x;
    default:
      return 0.0;
    }
  } 
  case 3 : { 
    switch (i) {
    case 0: return 0.5*SQR(0.5 - x);
    case 1: return 0.75 - SQR(x);
    case 2: return 0.5*SQR(0.5 + x);
    default:
      return 0.0;
    }
  case 4 : { 
    switch (i) {
    case 0: return ( 1.0+x*( -6.0+x*( 12.0-x* 8.0)))/48.0;
    case 1: return (23.0+x*(-30.0+x*(-12.0+x*24.0)))/48.0;
    case 2: return (23.0+x*( 30.0+x*(-12.0-x*24.0)))/48.0;
    case 3: return ( 1.0+x*(  6.0+x*( 12.0+x* 8.0)))/48.0;
    default:
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
      return 0.0;
    }
  }
  }}
  return 0.0;
}

__global__ void assign_charges(const double * const pos, const double * const q,
cufftDoubleComplex *mesh, const int m_size, const int cao, const double pos_shift, const
double hi) {
      /** id of the particle **/
      int id = blockIdx.x;
      /** position relative to the closest gird point **/
      double m_pos[3];
      /** index of the nearest mesh point **/
      int nmp_x, nmp_y, nmp_z;      

      m_pos[0] = pos[3*id + 0] * hi - pos_shift;
      m_pos[1] = pos[3*id + 1] * hi - pos_shift;
      m_pos[2] = pos[3*id + 2] * hi - pos_shift;

      nmp_x = (int) floor(m_pos[0] + 0.5);
      nmp_y = (int) floor(m_pos[1] + 0.5);
      nmp_z = (int) floor(m_pos[2] + 0.5);

      m_pos[0] -= nmp_x;
      m_pos[1] -= nmp_y;
      m_pos[2] -= nmp_z;

      nmp_x = wrap_index(nmp_x + threadIdx.x, m_size);
      nmp_y = wrap_index(nmp_y + threadIdx.y, m_size);
      nmp_z = wrap_index(nmp_z + threadIdx.z, m_size);

      /* printf("id %d, m { %d %d %d }: weight = %lf, nmp[] = (%d %d %d), pos[] = (%lf %lf %lf)\n", id, threadIdx.x, threadIdx.y, threadIdx.z, caf(threadIdx.x, m_pos[0], cao)*caf(threadIdx.y, m_pos[1], cao)*caf(threadIdx.z, m_pos[2], cao)*q[id], nmp_x, nmp_y, nmp_z, pos[0], pos[1], pos[2]); */

      atomicAdd( &(mesh[m_size*m_size*nmp_x +  m_size*nmp_y + nmp_z].x), caf(threadIdx.x, m_pos[0], cao)*caf(threadIdx.y, m_pos[1], cao)*caf(threadIdx.z, m_pos[2], cao)*q[id]);
}

__global__ void assign_forces(const double * const pos, const double * const q,
cufftDoubleComplex *mesh, const int m_size, const int cao, const double pos_shift, const
			      double hi, double *force, double prefactor) {
      /** id of the particle **/
      int id = blockIdx.x;
      /** position relative to the closest gird point **/
      double m_pos[3];
      /** index of the nearest mesh point **/
      int nmp_x, nmp_y, nmp_z;      

      m_pos[0] = pos[3*id + 0] * hi - pos_shift;
      m_pos[1] = pos[3*id + 1] * hi - pos_shift;
      m_pos[2] = pos[3*id + 2] * hi - pos_shift;

      nmp_x = (int) floor(m_pos[0] + 0.5);
      nmp_y = (int) floor(m_pos[1] + 0.5);
      nmp_z = (int) floor(m_pos[2] + 0.5);

      m_pos[0] -= nmp_x;
      m_pos[1] -= nmp_y;
      m_pos[2] -= nmp_z;

      nmp_x = wrap_index(nmp_x + threadIdx.x, m_size);
      nmp_y = wrap_index(nmp_y + threadIdx.y, m_size);
      nmp_z = wrap_index(nmp_z + threadIdx.z, m_size);

      /* printf("id %d, m { %d %d %d }: weight = %lf, nmp[] = (%d %d %d), pos[] = (%lf %lf %lf)\n", id, threadIdx.x, threadIdx.y, threadIdx.z, caf(threadIdx.x, m_pos[0], cao)*caf(threadIdx.y, m_pos[1], cao)*caf(threadIdx.z, m_pos[2], cao)*q[id], nmp_x, nmp_y, nmp_z, pos[0], pos[1], pos[2]); */

      atomicAdd( &(force[id]), prefactor*mesh[m_size*m_size*nmp_x +  m_size*nmp_y + nmp_z].x*caf(threadIdx.x, m_pos[0], cao)*caf(threadIdx.y, m_pos[1], cao)*caf(threadIdx.z, m_pos[2], cao)*q[id]);
}

__global__ void influence_function( double *G_hat, double box, int cao, int mesh, double alpha ) {
  int n[3];
  int linear_index;
  double nom[3] = { 0.0, 0.0, 0.0 }, dnom = 0.0;
  double fak = SQR( PI / alpha );
  int mx, my, mz;
  double box_i = 1./box;
  int nshift[3];
  int nmx, nmy, nmz;
  double zwi, nm2;
  double S1, S2, S3;

  n[0] = blockIdx.x;
  n[1] = blockIdx.y;
  n[2] = threadIdx.x;

  nshift[0] = (n[0] > mesh/2) ? n[0] - mesh : n[0];
  nshift[1] = (n[1] > mesh/2) ? n[1] - mesh : n[1];
  nshift[2] = (n[2] > mesh/2) ? n[2] - mesh : n[2];

  linear_index = SQR(mesh)*n[0] + mesh * n[1] + n[2];

  if( (n[0] == 0) && ( n[1] == 0) && n[2] == 0) {
    G_hat[linear_index] = 0.0;
    return;
  }

  if( (n[0] % (mesh/2) == 0)  && (n[1] % (mesh/2) == 0)  && (n[2] % (mesh/2) == 0)) {
    G_hat[linear_index] = 0.0;
    return;
  } 

  for ( mx = -P3M_BRILLOUIN; mx <= P3M_BRILLOUIN; mx++ ) {
    nmx = nshift[0] + mesh*mx;
    S1 = pow ( sinc ( box_i*nmx ), 2*cao );
    for ( my = -P3M_BRILLOUIN; my <= P3M_BRILLOUIN; my++ ) {
      nmy = nshift[1] + mesh*my;
      S2   = S1*pow ( sinc ( box_i*nmy ), 2*cao );
      for ( mz = -P3M_BRILLOUIN; mz <= P3M_BRILLOUIN; mz++ ) {
	nmz = nshift[2] + mesh*mz;
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

  printf("influence_function(%d %d %d) = %lf, nm2 = %lf, nm[] = (%d %d %d), nshift[] = (%d %d %d), dnom = %e\n",
	 n[0], n[1], n[2], zwi, nm2, nmx, nmy, nmz, nshift[0], nshift[1], nshift[2], dnom);
  
  G_hat[linear_index] = 2.0 * zwi / PI;

  return;
}

__global__ void apply_influence_function( cufftDoubleComplex *mesh, int mesh_size, double *G_hat ) {
  int linear_index = mesh_size*mesh_size*blockIdx.x + mesh_size * blockIdx.y + threadIdx.x;
  mesh[linear_index].x *= G_hat[linear_index];
  mesh[linear_index].y *= G_hat[linear_index];
}

__global__ void apply_diff_op( cufftDoubleComplex *mesh, const int mesh_size, cufftDoubleComplex *force_mesh,  const double box, const int dim ) {
  int linear_index = mesh_size*mesh_size*blockIdx.x + mesh_size * blockIdx.y + threadIdx.x;
  int n;

  switch( dim ) {
  case 0:
    n = blockIdx.x;
    break;
  case 1:
    n = blockIdx.y;
    break;
  case 2:
    n = threadIdx.x;
    break;
  }

  n = ( n == mesh_size/2 ) ? 0.0 : n;
  n = ( n > mesh_size/2) ? n - mesh_size : n;
 
  force_mesh[linear_index].x =  -2.0 * PI * n * mesh[linear_index].y / box;
  force_mesh[linear_index].y =   2.0 * PI * n * mesh[linear_index].x / box;
}

/* __global__ void assign_charges(const double * const pos, const double * const q, */
/* cufftDoubleComplex *mesh, const int m_size, const int cao, const double pos_shift, const */
/* double hi) { */

int main(void) {
  cufftHandle plan;
  cufftDoubleComplex *data, *force_mesh;

  const double pos_h[6] = {  4.0,  5.0,  5.0, 6.0, 5.0, 5.0 }, q_h[2] = { -1.0, 1.0 };
  double *pos_d, *q_d;
  const int cao = 7;
  const double box = 10.0;

  double forces_h[3][2], *forces_d;

  data_t *d;

  double *g_hat_d, *g_hat_h;
  cudaMalloc((void**)&g_hat_d, sizeof(double)*N*N*N);
  if (cudaGetLastError() != cudaSuccess){
    fprintf(stderr, "Cuda error: Failed to allocate\n");
    return 0;	
  }

  d = read_reference("test.dat");

  for(int i =0; i<d->n; ++i) {
    printf("pos %lf %lf %lf q %lf\n", d->pos[3*i+0], d->pos[3*i+1], d->pos[3*i+2], d->q[i]);
  }

  g_hat_h = (double *)malloc(N*N*N*sizeof(double));

  cufftDoubleComplex *data_h = (cufftDoubleComplex *) malloc( N * N * N * sizeof(cufftDoubleComplex));

  cudaMalloc((void**)&data, sizeof(cufftDoubleComplex)*N*N*N);
  if (cudaGetLastError() != cudaSuccess){
    fprintf(stderr, "Cuda error: Failed to allocate\n");
    return 0;	
  }

  cudaMalloc((void**)&force_mesh, sizeof(cufftDoubleComplex)*N*N*N);
  if (cudaGetLastError() != cudaSuccess){
    fprintf(stderr, "Cuda error: Failed to allocate\n");
    return 0;	
  }

  cudaMalloc((void**)&pos_d, sizeof(pos_h));
  cudaMalloc((void**)&q_d, sizeof(q_h));
  cudaMalloc((void**)&forces_d, sizeof(forces_h));

  cudaMemcpy( pos_d, pos_h, sizeof(pos_h), cudaMemcpyHostToDevice);
  cudaMemcpy( q_d, q_h, sizeof(q_h), cudaMemcpyHostToDevice);

  // prepare influence function
  dim3 blockDim(N, N, 1);
  dim3 thdDim( N, 1, 1);

  Influence_function_berechnen_ik( cao, N, box, 0.8, g_hat_h );

  cudaMemcpy( g_hat_d, g_hat_h, N*N*N*sizeof(double), cudaMemcpyHostToDevice);
\
  if (cudaThreadSynchronize() != cudaSuccess){
    fprintf(stderr, "Cuda error: Failed to synchronize\n");
    return 0;
  }

  cudaMemset( data, 0, N*N*N*sizeof(cufftDoubleComplex));
  
  dim3 caoBlock(cao, cao, cao);

  assign_charges<<<2, caoBlock>>>( pos_d, q_d, data, N, cao, (cao-1)/2.0, N/box);

  cudaThreadSynchronize();

  if (cufftPlan3d(&plan, N, N, N, CUFFT_Z2Z) != CUFFT_SUCCESS){
    fprintf(stderr, "CUFFT error: Plan creation failed");
    return 0;
  }

  if (cufftExecZ2Z(plan, data, data, CUFFT_FORWARD) != CUFFT_SUCCESS){
    fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed\n");
    return 0;
  }

  if (cudaThreadSynchronize() != cudaSuccess){
    fprintf(stderr, "Cuda error: Failed to synchronize\n");
    return 0;
  }

  apply_influence_function<<<blockDim, thdDim>>>( data, N, g_hat_d);

  for(int dim = 0; dim < 3; dim++) {
    if (cudaThreadSynchronize() != cudaSuccess){
      fprintf(stderr, "Cuda error: Failed to synchronize\n");
      return 0;
    }

    apply_diff_op<<<blockDim, thdDim>>>( data, N, force_mesh, box, dim);

    if (cudaThreadSynchronize() != cudaSuccess){
      fprintf(stderr, "Cuda error: Failed to synchronize diff_op\n");
      return 0;
    }

    /* Use the CUFFT plan to transform the signal in place. */
    if (cufftExecZ2Z(plan, force_mesh, force_mesh, CUFFT_INVERSE) != CUFFT_SUCCESS){
      fprintf(stderr, "CUFFT error: ExecZ2Z Backward failed\n");
      return 0;
    }

    if (cudaThreadSynchronize() != cudaSuccess){
      fprintf(stderr, "Cuda error: Failed to synchronize back\n");
      return 0;
    }

    cudaMemset(forces_d, 0, sizeof(forces_h));

    assign_forces<<< 2, caoBlock>>>( pos_d, q_d, force_mesh, N, cao, (cao-1)/2.0, N/box, forces_d, 1.0 / ( 2.0 *  box * box * box));

    cudaMemcpy( forces_h[dim], forces_d, sizeof(forces_h), cudaMemcpyDeviceToHost);
  }

  cudaMemcpy( data_h, data, sizeof(cufftDoubleComplex)*N*N*N, cudaMemcpyDeviceToHost);

  cufftDestroy(plan);
  cudaFree(data);

  printf("f_x = [%lf %lf]\n", forces_h[0][0], forces_h[0][1]); 
  printf("f_y = [%lf %lf]\n", forces_h[1][0], forces_h[1][1]); 
  printf("f_z = [%lf %lf]\n", forces_h[2][0], forces_h[2][1]); 

}

