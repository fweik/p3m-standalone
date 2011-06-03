#ifndef P3M_H
#define P3M_H

#define MaxInterpol (2*100096)
#define Maxip 6

/* Pi, weil man's so oft braucht: */
#define PI 3.14159265358979323846264

#define SQR(A) ((A)*(A))

//Pointer to Force Arrays for kspace
double *Fx_K, *Fy_K, *Fz_K;
//Pointer to particle charges
double *Q;
//Pointer to particle possitions;
double *xS, *yS, *zS;
//Pointer to influence function
double *G_hat;
// Mesh size
int Mesh;
// Box length
double Len, Leni;
// charge assignment order - 1
int ip;
// Pointer to mesh charge density
double *Qmesh;
// differentials
double *dQdx, *dQdy, *dQdz;
// Pointer to differential operator
double *Dn;
//
/* Speichert die Interpolation des Ladungszuordnungspolynoms: */
double LadInt[Maxip+1][2*MaxInterpol+1];
/* Speichert die Interpolation der Ableitung des Ladungszuordnungspolynoms: */
double LadInt_[Maxip+1][2*MaxInterpol+1];

double Q2;

double E_Coulomb_Dipol;
double E_Coulomb_Self;
double E_Coulomb_Impuls_Summe;
double E_Coulomb_Real_Summe;

// Internal variables
int *Gx, *Gy, *Gz;
double *nshift;

#endif
