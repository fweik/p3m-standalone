#ifndef P3M_H
#define P3M_H

#define MaxInterpol (2*100096)
#define Maxip 6

// #define SINGLE_PREC
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

/* Pi, weil man's so oft braucht: */
#define PI 3.14159265358979323846264

#define SQR(A) ((A)*(A))

//Pointer to Force Arrays for kspace
FLOAT_TYPE *Fx_K, *Fy_K, *Fz_K;
//Pointer to particle charges
FLOAT_TYPE *Q;
//Pointer to particle possitions;
FLOAT_TYPE *xS, *yS, *zS;
//Pointer to influence function
FLOAT_TYPE *G_hat;
// Mesh size
int Mesh;
// Box length
FLOAT_TYPE Len, Leni;
// charge assignment order - 1
int ip;
// charge assignment order
int cao, cao3;
// Pointer to mesh charge density
FLOAT_TYPE *Qmesh;
// differentials
FLOAT_TYPE *dQdx, *dQdy, *dQdz;
// Pointer to differential operator
FLOAT_TYPE *Dn;
//
/* Speichert die Interpolation des Ladungszuordnungspolynoms: */
FLOAT_TYPE LadInt[Maxip+1][2*MaxInterpol+1];
/* Speichert die Interpolation der Ableitung des Ladungszuordnungspolynoms: */
FLOAT_TYPE LadInt_[Maxip+1][2*MaxInterpol+1];

FLOAT_TYPE Q2;

FLOAT_TYPE E_Coulomb_Dipol;
FLOAT_TYPE E_Coulomb_Self;
FLOAT_TYPE E_Coulomb_Impuls_Summe;
FLOAT_TYPE E_Coulomb_Real_Summe;

// Internal variables
int *Gx, *Gy, *Gz;
FLOAT_TYPE *nshift;

// Charge assignment

int *ca_ind[2];
FLOAT_TYPE *cf[2];

#endif
