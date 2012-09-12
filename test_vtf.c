#include "types.h"
#include "io.h"
#include "generate_system.h"

int main(void) {
  system_t *s = generate_system( FORM_FACTOR_SEPERATED_DIPOLE, 30, 10.0, 1.0);

  write_vtf( "test.vtf", s);

  return 0;
}
