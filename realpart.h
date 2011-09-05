#ifndef REALPART_H
#define REALPART_H

#include "types.h"

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

// Error of the realspace part

FLOAT_TYPE Realspace_error( const system_t *, const parameters_t * );

#endif
