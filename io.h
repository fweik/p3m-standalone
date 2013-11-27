#ifndef IO_H
#define IO_H

#include "types.h"

void Read_exact_forces(system_t *s, char *);
system_t *Read_system(parameters_t *, char *);
void Write_exact_forces(system_t *, char *);
void Write_system(system_t *, char *);
void Write_system_cuda( system_t *s, parameters_t *p, char *filename);

void write_vtf(char *filename, system_t *s);

#endif
