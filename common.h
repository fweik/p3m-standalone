#pragma once

#ifndef COMMON_H
#define COMMON_H

#include <stdlib.h>

// Set floating point precision

#define DOUBLE_PREC

#ifdef SINGLE_PREC
  #define FLOAT_FORMAT "%.8f"
  #define FLOAT_TYPE float
#endif

#ifdef DOUBLE_PREC 
  #define FLOAT_FORMAT "%.15f"
  #define FLOAT_TYPE double
#endif

#ifdef QUAD_PREC
  #define FLOAT_TYPE __float128
  #define FLOAT_FORMAT "%.35q"
#endif

// !3
#define PI 3.14159265358979323846264

#define SQR(A) ((A)*(A))

// Method ids

enum {
  METHOD_P3M_ik = 0,
  METHOD_P3M_ik_i = 1,
  METHOD_P3M_ad = 2,
  METHOD_P3M_ad_i = 3,
  METHOD_EWALD = 4
};

// Container type for arrays of 3d-vectors
// each component holds a pointer to an array
// of the values for that direction.
// Fields contains the same pointsers as array
// for convinience.

typedef struct {
  FLOAT_TYPE *x;
  FLOAT_TYPE *y;
  FLOAT_TYPE *z;
  FLOAT_TYPE **fields;
  } vector_array_t;

// Struct to hold the forces
  
typedef struct {
  // Total forces
  vector_array_t f;
  // Kspace part of force
  vector_array_t f_k;
  // Realspace part of force
  vector_array_t f_r;
} forces_t;
  
// Struct to hold the reference forces

typedef struct {
  vector_array_t f;
  vector_array_t f_k;
} reference_forces_t;

// Struct to hold general system and particle data

typedef
  struct {
    // box length
    FLOAT_TYPE length;
    // number of particles
    int        nparticles;
    // particle positions;
    vector_array_t p; 
    // particle forces
    vector_array_t f;
    // longrange part of forces
    vector_array_t f_k;
    // shortrange part of forces
    vector_array_t f_r;
    // charges of the particles
    FLOAT_TYPE *q;
    // sum of the squares of the particle charges
    FLOAT_TYPE q2;
    // the reference forces
    reference_forces_t reference;
    // dielectric constat of the environment at infinity
    FLOAT_TYPE epsilon;
} system_t;

void Init_system(system_t *);
void *Init_array(int, size_t);
void Init_vector_array(vector_array_t *, int);
void Free_vector_array(vector_array_t *);
#endif
