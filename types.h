#ifndef TYPES_H
#define TYPES_H

#include <fftw3.h>

// Set floating point precision

#ifdef DETAILED_TIMINGS
extern int __detailed_timings;
#endif

// Limit valgrind profiling to interesting part
// #define __VALGRIND_PROFILE_KSPACE_ONLY

#define DOUBLE_PREC 
//#define LONG_DOUBLE_PREC

#ifdef SINGLE_PREC
#define DIGITS 9
#define FLOAT_TYPE float
#define FLOAT_CAST (double)
#define EXP expf
#define SIN sinf
#define COS cosf
#define SQRT sqrtf
#define ERFC erfcf
#define FLOAT_ABS fabsf
#define FFTW_FREE fftwf_free
#define FFTW_MALLOC fftwf_malloc
#define FFTW_EXECUTE fftwf_execute
#define FFTW_COMPLEX fftwf_complex
#define FFTW_PLAN_DFT_3D fftwf_plan_dft_3d
#define FFTW_PLAN fftwf_plan
#define FFTW_DESTROY_PLAN fftwf_destroy_plan
#define ROUND roundf
#define FLOOR floorf
#endif

#ifdef DOUBLE_PREC
#define DIGITS 17
#define FLOAT_TYPE double
#define FLOAT_CAST
#define EXP exp
#define SIN sin
#define COS cos
#define SQRT sqrt
#define ERFC erfc
#define FLOAT_ABS fabs
#define FFTW_FREE fftw_free
#define FFTW_MALLOC fftw_malloc
#define FFTW_EXECUTE fftw_execute
#define FFTW_COMPLEX fftw_complex
#define FFTW_PLAN_DFT_3D fftw_plan_dft_3d
#define FFTW_PLAN_DFT_R2C_3D fftw_plan_dft_r2c_3d
#define FFTW_PLAN_DFT_C2R_3D fftw_plan_dft_c2r_3d
#define FFTW_PLAN fftw_plan
#define FFTW_DESTROY_PLAN fftw_destroy_plan
#define ROUND round
#define FLOOR floor
#define LOG log
#endif

#ifdef LONG_DOUBLE_PREC
#define DIGITS 22
#define FLOAT_TYPE long double
#define FLOAT_FORMAT "%.35q"
#define FLOAT_CAST (double)
#define EXP expl
#define SIN sinl
#define COS cosl
#define SQRT sqrtl
#define ERFC erfcl
#define FLOAT_ABS fabsl
#define FFTW_FREE fftwl_free
#define FFTW_MALLOC fftwl_malloc
#define FFTW_EXECUTE fftwl_execute
#define FFTW_COMPLEX fftwl_complex
#define FFTW_PLAN_DFT_3D fftwl_plan_dft_3d
#define FFTW_PLAN_DFT_R2C_3D fftw_plan_dft_r2c_3d
#define FFTW_PLAN_DFT_C2R_3D fftw_plan_dft_c2r_3d
#define FFTW_PLAN fftwl_plan
#define FFTW_DESTROY_PLAN fftwl_destroy_plan
#define ROUND roundl
#define FLOOR floorl
#define LOG logl
#endif

// !3
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062

#define BLIST_STEP 100

#define SQR(A) ((A)*(A))

// Method ids

enum {
    METHOD_P3M_ik = 0,
    METHOD_P3M_ik_i = 1,
    METHOD_P3M_ad = 2,
    METHOD_P3M_ad_i = 3,
    METHOD_EWALD = 4,
    METHOD_P3M_ik_cuda = 5,
    METHOD_P3M_ik_r = 6
};

// Container type for arrays of 3d-vectors
// each component holds a pointer to an array
// of the values for that direction.
// Fields contains the same pointers as array
// for convinience.

typedef struct {
  int size;
  FLOAT_TYPE *x;
  FLOAT_TYPE *y;
  FLOAT_TYPE *z;
  FLOAT_TYPE **fields;
} vector_array_t;


typedef struct {
  size_t size;
  size_t bufsize;
  void *data;
} buffered_list_t;

typedef struct {
  int size;
  FLOAT_TYPE *x;
  FLOAT_TYPE *y;
  FLOAT_TYPE *z;
  buffered_list_t **data;
  FLOAT_TYPE **fields;
} bvector_array_t;

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
    // particle veloceties
    vector_array_t *v;
    // the reference forces
    forces_t *reference;
    // dielectric constat of the environment at infinity
    FLOAT_TYPE epsilon;
    // energy of the system
    FLOAT_TYPE energy;
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

typedef struct {
  // cao order of the interpolation
  int cao;
  // Actual interpolation data
  FLOAT_TYPE **interpol;
  // Interpolation of the derivative
  FLOAT_TYPE **interpol_d;
  // array function pointers to the FT of the CA-function
  FLOAT_TYPE (*U_hat)(int, FLOAT_TYPE);
} interpolation_t;

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
    // Struct for interpolated charge assignment function
    interpolation_t *inter;
    // fftw plans
    //number of plans
    int forward_plans;
    int backward_plans;
    // actual plans
    FFTW_PLAN forward_plan[3];
    FFTW_PLAN backward_plan[3];
    // neighbor list for real space calculation
    neighbor_list_t *neighbor_list;
    // Self forces corrections
    FLOAT_TYPE *self_force_corrections;
    void *method_data;
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
    METHOD_FLAG_self_force_correction = 128, // Method need self force correction
};

// Common flags for all p3m methods for convinience
#define METHOD_FLAG_P3M (METHOD_FLAG_nshift | METHOD_FLAG_G_hat | METHOD_FLAG_Qmesh | METHOD_FLAG_ca )
// methode type

typedef struct {
    int  method_id;
    const char *method_name;
    const char *method_name_short;
    char flags;
    data_t * ( *Init ) ( system_t *, parameters_t * );
    void ( *Influence_function ) ( system_t *, parameters_t *, data_t * );
    void ( *Kspace_force ) ( system_t *, parameters_t *, data_t *, forces_t * );
    FLOAT_TYPE ( *Error ) ( system_t *, parameters_t * );
    FLOAT_TYPE ( *Error_k ) ( system_t *, parameters_t * );
} method_t;

// Function pointer types

typedef FLOAT_TYPE (*R3_to_R)(int, int, int, system_t *s, parameters_t *p);
typedef FLOAT_TYPE (*R_to_R)(FLOAT_TYPE);

#endif
