#pragma once

#ifndef DIPOL_H
#define DIPOL_H

/* Calculate the dipol correction to the Ewald
 * sum in the case that the dielectric constant
 * of the surronding environment is < \infty
 */

#include "common.h"

void Dipol(system_t *s);

#endif