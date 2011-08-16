#ifndef P3M_H
#define P3M_H

#include "common.h"
#include "p3m-common.h"

#define MaxInterpol (2*100096)
#define Maxip 6

/* Speichert die Interpolation des Ladungszuordnungspolynoms: */
FLOAT_TYPE LadInt[Maxip+1][2*MaxInterpol+1];
/* Speichert die Interpolation der Ableitung des Ladungszuordnungspolynoms: */
FLOAT_TYPE LadInt_[Maxip+1][2*MaxInterpol+1];

// Charge assignment

int *ca_ind[2];
FLOAT_TYPE *cf[2];



// Struct holding p3m parameters.

typedef struct {
  FLOAT_TYPE alpha;
  FLOAT_TYPE rcut;
  FLOAT_TYPE prefactor;
  int        mesh;
  int        ip;
  int        cao;
  int        cao3;
} p3m_parameters_t;

// Struct holding p3m data.

typedef struct {
  int mesh;
  FLOAT_TYPE *G_hat;
  FLOAT_TYPE *Qmesh;
  vector_array_t Fmesh;
  FLOAT_TYPE *nshift;
  FLOAT_TYPE *Dn;
  FLOAT_TYPE *dQdx[2], *dQdy[2], *dQdz[2];
} p3m_data_t;

// Flags for method_t

enum {
  P3M_FLAG_ik = 1,
  P3M_FLAG_ad = 2,
  P3M_FLAG_interlaced = 4,
  P3M_FLAG_nshift = 8
};

// methode type

typedef struct {
  int  method_id;
  const char *method_name;
  char flags;
  void (*Init)(system_t *, p3m_parameters_t *);
  void (*Influence_function)(system_t *, p3m_parameters_t *);
  void (*Kspace_force)(system_t *, p3m_parameters_t *);
  double (*Error)(system_t *, p3m_parameters_t *);
} method_t;

#endif
