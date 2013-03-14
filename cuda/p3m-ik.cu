#include <cufft.h>

#define PI 3.14159265359

#include "p3m-ik-cuda.h"

#define SQR(A) ((A)*(A))
#define P3M_BRILLOUIN 1

using namespace std;

typedef struct {
  cufftHandle plan;
  cufftDoubleComplex *charge_mesh;
  cufftDoubleComplex *force_mesh;
  double *g_hat_d;
  double *pos_d;
  double *q_d;
  double *forces_d;
  p3m_cuda_data_t *d;
} p3m_cuda_state_t;

void p3m_ik_cuda_free( p3m_cuda_data_t *d) {
  p3m_cuda_state_t *p3m_cuda_state = (p3m_cuda_state_t *) d->s;  

  cufftDestroy(p3m_cuda_state->plan);
 
  cudaFree(p3m_cuda_state->charge_mesh);
  cudaFree(p3m_cuda_state->force_mesh);

  cudaFree(p3m_cuda_state->g_hat_d);
  cudaFree(p3m_cuda_state->q_d);
  cudaFree(p3m_cuda_state->pos_d);
  cudaFree(p3m_cuda_state->forces_d);
}

__device__ __host__ inline static double sinc(double d)
{
  double PId = PI*d;
  return (d == 0.0) ? 1.0 : sin(PId)/PId;
}

void static Aliasing_sums_ik ( int cao, double box, double alpha, int mesh, int NX, int NY, int NZ,
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
void static Influence_function_berechnen_ik ( int cao, int mesh, double box, double alpha, double *G_hat ) {

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

__device__ double caf(int i, double x, int cao_value) {
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

      /* printf("id %d, m { %d %d %d }: weight = %lf, nmp[] = (%d %d %d), pos[] = (%lf %lf %lf)\n", id, threadIdx.x, threadIdx.y, threadIdx.z, caf(threadIdx.x, m_pos[0], cao)*caf(threadIdx.y, m_pos[1], cao)*caf(threadIdx.z, m_pos[2], cao)*q[id], nmp_x, nmp_y, nmp_z, m_pos[0], m_pos[1], m_pos[2]); */

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

      atomicAdd( &(force[id]), -prefactor*mesh[m_size*m_size*nmp_x +  m_size*nmp_y + nmp_z].x*caf(threadIdx.x, m_pos[0], cao)*caf(threadIdx.y, m_pos[1], cao)*caf(threadIdx.z, m_pos[2], cao)*q[id]);
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
 
void p3m_ik_cuda_init( p3m_cuda_data_t *d ) {
  puts("p3m_ik_cuda_init():");
  double *g_hat_h = (double *)malloc(d->mesh*d->mesh*d->mesh*sizeof(double));
  p3m_cuda_state_t *p3m_cuda_state = (p3m_cuda_state_t *) malloc(sizeof(p3m_cuda_state_t));

  puts("Allocating g_hat_d.");
  cudaMalloc((void**)&(p3m_cuda_state->g_hat_d), sizeof(double)*d->mesh*d->mesh*d->mesh);
  if (cudaGetLastError() != cudaSuccess){
    fprintf(stderr, "p3m_cuda: Failed to allocate\n");
  }

  puts("Allocating charge_mesh.");
  cudaMalloc((void**)&(p3m_cuda_state->charge_mesh), sizeof(cufftDoubleComplex)*d->mesh*d->mesh*d->mesh);
  if (cudaGetLastError() != cudaSuccess){
    fprintf(stderr, "p3m_cuda: Failed to allocate\n");
  }

  puts("Allocating foce_mesh.");
  cudaMalloc((void**)&(p3m_cuda_state->force_mesh), sizeof(cufftDoubleComplex)*d->mesh*d->mesh*d->mesh);
  if (cudaGetLastError() != cudaSuccess){
    fprintf(stderr, "p3m_cuda: Failed to allocate\n");
  }


  cudaMalloc((void**)&(p3m_cuda_state->pos_d), 3*d->n*sizeof(double));
  if (cudaGetLastError() != cudaSuccess){
    fprintf(stderr, "p3m_cuda: Failed to allocate\n");
  }
  cudaMalloc((void**)&(p3m_cuda_state->q_d), d->n*sizeof(double));
  if (cudaGetLastError() != cudaSuccess){
    fprintf(stderr, "p3m_cuda: Failed to allocate\n");
  }
  cudaMalloc((void**)&(p3m_cuda_state->forces_d), d->n*sizeof(double));
  if (cudaGetLastError() != cudaSuccess){
    fprintf(stderr, "p3m_cuda: Failed to allocate\n");
  }

  Influence_function_berechnen_ik( d->cao, d->mesh, d->box, d->alpha, g_hat_h );  

  cudaMemcpy( p3m_cuda_state->g_hat_d, g_hat_h, d->mesh*d->mesh*d->mesh*sizeof(double), cudaMemcpyHostToDevice);

  if (cufftPlan3d(&(p3m_cuda_state->plan), d->mesh, d->mesh, d->mesh, CUFFT_Z2Z) != CUFFT_SUCCESS){
    fprintf(stderr, "CUFFT error: Plan creation failed");
  }
}

void p3m_ik_cuda(p3m_cuda_data_t *d) {
  p3m_cuda_state_t *p3m_cuda_state = (p3m_cuda_state_t *) d->s;

  cudaMemcpy( p3m_cuda_state->pos_d, d->pos, 3*d->n*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy( p3m_cuda_state->q_d, d->q, d->n*sizeof(double), cudaMemcpyHostToDevice);

  dim3 blockDim(d->mesh, d->mesh, 1);
  dim3 thdDim( d->mesh, 1, 1);

  if (cudaThreadSynchronize() != cudaSuccess){
    fprintf(stderr, "Cuda error: Failed to synchronize\n");
    return;
  }

  cudaMemset( p3m_cuda_state->charge_mesh, 0, d->mesh*d->mesh*d->mesh*sizeof(cufftDoubleComplex));
  
  dim3 caoBlock(d->cao, d->cao, d->cao);

  assign_charges<<<d->n, caoBlock>>>( p3m_cuda_state->pos_d, p3m_cuda_state->q_d, p3m_cuda_state->charge_mesh, d->mesh, d->cao,(double)((d->cao-1)/2), d->mesh/d->box);

  cudaThreadSynchronize();

  if (cufftExecZ2Z(p3m_cuda_state->plan, p3m_cuda_state->charge_mesh, p3m_cuda_state->charge_mesh, CUFFT_FORWARD) != CUFFT_SUCCESS){
    fprintf(stderr, "CUFFT error: ExecZ2Z Forward failed\n");
    return;
  }

  if (cudaThreadSynchronize() != cudaSuccess){
    fprintf(stderr, "Cuda error: Failed to synchronize\n");
    return;
  }

  apply_influence_function<<<blockDim, thdDim>>>( p3m_cuda_state->charge_mesh, d->mesh, p3m_cuda_state->g_hat_d);

  for(int dim = 0; dim < 3; dim++) {
    if (cudaThreadSynchronize() != cudaSuccess){
      fprintf(stderr, "Cuda error: Failed to synchronize\n");
      return;
    }

    apply_diff_op<<<blockDim, thdDim>>>( p3m_cuda_state->charge_mesh, d->mesh, p3m_cuda_state->force_mesh, d->box, dim);

    if (cudaThreadSynchronize() != cudaSuccess){
      fprintf(stderr, "Cuda error: Failed to synchronize diff_op\n");
      return;
    }

    /* Use the CUFFT plan to transform the signal in place. */
    if (cufftExecZ2Z(p3m_cuda_state->plan, p3m_cuda_state->force_mesh, p3m_cuda_state->force_mesh, CUFFT_INVERSE) != CUFFT_SUCCESS){
      fprintf(stderr, "CUFFT error: ExecZ2Z Backward failed\n");
      return;
    }

    if (cudaThreadSynchronize() != cudaSuccess){
      fprintf(stderr, "Cuda error: Failed to synchronize back\n");
      return;
    }

    cudaMemset(p3m_cuda_state->forces_d, 0, d->n*sizeof(double));

    assign_forces<<< d->n, caoBlock>>>( p3m_cuda_state->pos_d, p3m_cuda_state->q_d, p3m_cuda_state->force_mesh, d->mesh, d->cao, (double)((d->cao-1)/2), d->mesh/d->box, p3m_cuda_state->forces_d, 1.0 / ( 2.0 *  d->box * d->box * d->box));

    cudaMemcpy( d->f[dim], p3m_cuda_state->forces_d, d->n*sizeof(double), cudaMemcpyDeviceToHost);
  }
}

