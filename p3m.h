#pragma once

#ifndef P3M_H
#define P3M_H

#include "common.h"

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
    P3M_FLAG_none = 0,
    P3M_FLAG_ik = 1, // Method uses ik-diff
    P3M_FLAG_ad = 2, // Method uses ad
    P3M_FLAG_interlaced = 4, // Method is interlaced
    P3M_FLAG_nshift = 8, // Method needs precalculated shifted k-values
    P3M_FLAG_G_hat = 16, // Method uses influence function
    P3M_FLAG_Qmesh = 32, // Method needs charge mesh
    P3M_FLAG_ca = 64 // Method uses charge assignment
};

// methode type

typedef struct {
    int  method_id;
    const char *method_name;
    char flags;
    p3m_data_t * ( *Init ) ( system_t *, p3m_parameters_t * );
    void ( *Influence_function ) ( system_t *, p3m_parameters_t *, p3m_data_t * );
    void ( *Kspace_force ) ( system_t *, p3m_parameters_t *, p3m_data_t * );
    double ( *Error ) ( system_t *, p3m_parameters_t * );
} method_t;

#endif
