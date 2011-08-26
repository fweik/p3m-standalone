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

// Struct holding p3m parameters.

typedef struct {
    FLOAT_TYPE alpha;
    FLOAT_TYPE rcut;
    FLOAT_TYPE prefactor;
    int        mesh;
    int        ip;
    int        cao;
    int        cao3;
} parameters_t;

// Struct holding p3m data.

typedef struct {
    // Mesh size the struct is initialized for
    int mesh;
    // Influence function
    FLOAT_TYPE *G_hat;
    // Charge mesh
    FLOAT_TYPE *Qmesh;
    // Force mesh for k space differentiation
    vector_array_t Fmesh;
    // Shifted kvectors (fftw convention)
    FLOAT_TYPE *nshift;
    // Fourier coefficients of the differential operator
    FLOAT_TYPE *Dn;
    // Derivatives of the charge assignment function for analytical differentiation
    FLOAT_TYPE *dQdx[2], *dQdy[2], *dQdz[2];
    // Cache for charge assignment
    int *ca_ind[2];
    FLOAT_TYPE *cf[2];
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
    METHOD_FLAG_ca = 64 // Method uses charge assignment
};

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
