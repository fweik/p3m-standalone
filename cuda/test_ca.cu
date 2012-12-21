#include "assignment.h"

#include <stdlib.h>

#include <gsl/gsl_rng.h>

#include <cufft.h>

#define MESH_SIZE 32
#define PARTICLES 8
#define CAO 5
#define BOX 15.0

#define FLOAT_TYPE float



void generate_random_system(int size, FLOAT_TYPE box, FLOAT_TYPE max_charge, FLOAT_TYPE *pos, FLOAT_TYPE *q) {
  int i,j;
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

  for(i=0;i<size;i++) {
    for(j=0;j<3;j++) {
      pos[3*i+j] = box*gsl_rng_uniform(rng);
    }
    q[i] = (1.0 - 2.0 * (i%2)) * max_charge;
  }

  gsl_rng_free (rng);
}


int main(void) {

  cufftComplex *h_mesh = (cufftComplex *)malloc(MESH_SIZE*MESH_SIZE*MESH_SIZE*sizeof(cufftComplex));
  FLOAT_TYPE *h_part = (FLOAT_TYPE *)malloc(3*PARTICLES*sizeof(FLOAT_TYPE));
  FLOAT_TYPE *h_q    = (FLOAT_TYPE *)malloc(PARTICLES*sizeof(FLOAT_TYPE));

  FLOAT_TYPE *d_mesh, *d_part, *d_q;

  FLOAT_TYPE h = BOX / MESH_SIZE;

  FLOAT_TYPE sum = 0.0;

  cufftHandle plan;

  generate_random_system( PARTICLES, BOX, 1.0, h_part, h_q);

  for(int i = 0; i< PARTICLES; i++)
    printf("part %d, q = %f, pos = (%f %f %f)\n", i, h_q[i], h_part[3*i+0], h_part[3*i+1], h_part[3*i+2]);

  
  cudaMalloc((void**)&d_mesh, MESH_SIZE*MESH_SIZE*MESH_SIZE*sizeof(cufftComplex));
  cudaMalloc((void**)&d_part, 3*PARTICLES*sizeof(FLOAT_TYPE));
  cudaMalloc((void**)&d_q, PARTICLES*sizeof(FLOAT_TYPE));

  cudaMemcpy( d_part, h_part, 3*PARTICLES*sizeof(FLOAT_TYPE), cudaMemcpyHostToDevice);
  cudaMemcpy( d_q, h_q, PARTICLES*sizeof(FLOAT_TYPE), cudaMemcpyHostToDevice);
  cudaMemset( d_mesh, 0, MESH_SIZE*MESH_SIZE*MESH_SIZE*sizeof(cufftComplex));

  dim3 grid(PARTICLES,1), block(CAO,CAO,CAO);

  assign_charges<<<grid, block>>>( d_part, d_q, d_mesh, MESH_SIZE, CAO, (FLOAT_TYPE) ((CAO-1)/2), 1./h);

  cudaMemcpy( h_mesh, d_mesh, MESH_SIZE*MESH_SIZE*MESH_SIZE*sizeof(cufftComplex), cudaMemcpyDeviceToHost);

  for(int i = 0; i< MESH_SIZE*MESH_SIZE*MESH_SIZE; i++) {
    printf("%d %f\n", i, h_mesh[i].x);
    sum += h_mesh[i].x;
  }
  printf("sum %f\n", sum);
}

