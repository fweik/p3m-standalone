#ifndef IO_H
#define IO_H

#include "types.h"

void Read_exact_forces(system_t *s, char *);
system_t *Read_system(parameters_t *, char *);
void Write_exact_forces(system_t *, char *);

#endif
