#include "io.h"
#include "p3m-common.h"
#include "common.h"
#include "generate_system.h"

int main(void) {
  double box = 10.;
  int particles = 1000;
  system_t *system = generate_system(SYSTEM_RANDOM, particles, box, 1. );
  parameters_t parameters;
  double acc;

  acc = Calculate_reference_forces( system, &parameters );
  printf("accuracy %e\n", acc);
}

