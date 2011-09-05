#ifndef TYPES_H
#define TYPES_H

#include <fftw3.h>

// Set floating point precision

#ifdef DETAILED_TIMINGS
extern int __detailed_timings;
#endif

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
    METHOD_EWALD = 4,
    METHOD_GREENS_ik = 5
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
    vector_array_t *f;
    // Kspace part of force
    vector_array_t *f_k;
    // Realspace part of force
    vector_array_t *f_r;
} forces_t;

// Struct to hold general system and particle data

typedef
struct {
    // box length
    FLOAT_TYPE length;
    // number of particles
    int        nparticles;
    // particle positions;
    vector_array_t *p;
    // charges of the particles
    FLOAT_TYPE *q;
    // sum of the squares of the particle charges
    FLOAT_TYPE q2;
    // the reference forces
    forces_t *reference;
    // dielectric constat of the environment at infinity
    FLOAT_TYPE epsilon;
} system_t;

// struct for the neighbor list to speed up the real part calculation

typedef 
struct {
  // number of particles in list
  int n;
  // postitions of neighbors
  vector_array_t *p;
  // charges of neighbors
  FLOAT_TYPE *q;
  // ids in system array of neighbors
  int *id;
} neighbor_list_t;


// Struct holding method parameters.

typedef struct {
    FLOAT_TYPE alpha;
    FLOAT_TYPE rcut;
    FLOAT_TYPE prefactor;
    int        mesh;
    int        ip;
    int        cao;
    int        cao3;
    FLOAT_TYPE precision;
} parameters_t;

// Struct holding method data.

typedef struct {
    // Mesh size the struct is initialized for
    int mesh;
    // Influence function
    FLOAT_TYPE *G_hat;
    // Charge mesh
    FLOAT_TYPE *Qmesh;
    // Force mesh for k space differentiation
    vector_array_t *Fmesh;
    // Shifted kvectors (fftw convention)
    FLOAT_TYPE *nshift;
    // Fourier coefficients of the differential operator
    FLOAT_TYPE *Dn;
    // Derivatives of the charge assignment function for analytical differentiation
    FLOAT_TYPE *dQdx[2], *dQdy[2], *dQdz[2];
    // Cache for charge assignment
    int *ca_ind[2];
    FLOAT_TYPE *cf[2];
    // Array for interpolated charge assignment function
    FLOAT_TYPE **LadInt;
    FLOAT_TYPE **LadInt_;
    // fftw plans
    //number of plans
    int forward_plans;
    int backward_plans;
    // actual plans
    fftw_plan forward_plan[3];
    fftw_plan backward_plan[3];
    // neighbor list for real space calculation
    neighbor_list_t *neighbor_list;
} data_t;

// Flags for method_t

enum {
    METHOD_FLAG_none = 0,
    METHOD_FLAG_ik = 1, // Method uses ik-diff
    METHOD_FLAG_ad = 2, // Method uses ad
    METHOD_FLAG_interlaced = 4, // Method is interlaced
    METHOD_FLAG_nshift = 8, // Method needs precalculated shifted k-values
    METHOD_FLAG_G_hat = 16, // Method uses influence function
    METHOD_FLAG_Qmesh = 32, // Method needs charge mesh
    METHOD_FLAG_ca = 64, // Method uses charge assignment
};

// Common flags for all p3m methods for convinience
#define METHOD_FLAG_P3M (METHOD_FLAG_nshift | METHOD_FLAG_G_hat | METHOD_FLAG_Qmesh | METHOD_FLAG_ca )
// methode type

typedef struct {
    int  method_id;
    const char *method_name;
    char flags;
    data_t * ( *Init ) ( system_t *, parameters_t * );
    void ( *Influence_function ) ( system_t *, parameters_t *, data_t * );
    void ( *Kspace_force ) ( system_t *, parameters_t *, data_t *, forces_t * );
    double ( *Error ) ( system_t *, parameters_t * );
} method_t;


#endif
