#include "assignment.h"

#include <stdlib.h>

#include <gsl/gsl_rng.h>

#define mesh_size 32
#define particles 512
#define cao 5
#define box 15.0

int main(void) {
  float *h_mesh = (float *)malloc(mesh_size*sizeof(float));
  float *h_part = (float *)malloc(3*particles*sizeof(float));
  float *h_q    = (float *)malloc(particles*sizeof(float));

  float *d_mesh, *d_part, *d_q;

  float h = box / mesh_size;

  for(int i = 0; i<particles; i++) {

  }

  cudaMalloc((void**)&d_mesh, mesh_size*sizeof(float));
  cudaMalloc((void**)&d_part, 3*particles*sizeof(float));
  cudaMalloc((void**)&d_q, particles**sizeof(float));
}

