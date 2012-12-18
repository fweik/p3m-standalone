#ifndef REALPART_H
#define REALPART_H

#include "types.h"
#include "domain-decomposition.h"

/* Calculate Real space part of the Ewald sum.
 * Used in all methods
 */


// functions for neighbor list algorithm

// Build particle neighbor list. WARNING: This is O(n^2).
void Init_neighborlist(system_t *, parameters_t *, data_t * );
// Calculate realpart of forces using list.
void Realpart_neighborlist(system_t *, parameters_t *, data_t *, forces_t *);
// free neighor list
void Free_neighborlist(data_t *);

// functions for n2 algorithm

void Realteil(system_t *, parameters_t *, forces_t *);

FLOAT_TYPE Realpart_corr_error(FLOAT_TYPE rcut, FLOAT_TYPE alpha);

// Error of the realspace part

FLOAT_TYPE Realspace_error( const system_t *, const parameters_t * );

// Count neighbor pairs in s
int *count_neighbors( system_t *s, parameters_t *p );

#endif
