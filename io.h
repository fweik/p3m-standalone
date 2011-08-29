#ifndef IO_H
#define IO_H

#include "p3m.h"

void Exakte_Werte_einlesen(system_t *s, char *);
system_t *Daten_einlesen(parameters_t *, char *);
void Write_exact_forces(system_t *, char *);

#endif
