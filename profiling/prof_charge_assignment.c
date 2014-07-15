#include <stdio.h>
#include <mpi.h>

#include "../types.h"
#include "../charge-assign.h"
#include "../generate_system.h"
#include "../p3m-common.h"
#include "../p3m-ad.h"
#include "../p3m-ik.h"

int main(int argc, char **argv) {
  system_t *s = generate_system(0, atoi(argv[1]), 10.0, 1.0);
  parameters_t p;
  FLOAT_TYPE t, t_ad = 0.0, t_ik = 0.0;
  data_t *d_ad;
  data_t *d_ik;
  forces_t *f;

  f = Init_forces(s->nparticles);

  p.mesh = atoi(argv[2]);
  p.cao = atoi(argv[3]);
  p.ip = p.cao-1;
  p.cao3 = p.cao * SQR(p.cao); 
  p.rcut = 3.0; p.alpha = 0.5;

  d_ik = Init_ik( s, &p );

  t = MPI_Wtime();
  assign_charge( s, &p, d_ik, 0 );
  t -= MPI_Wtime();

  t_ik += t;

  printf("Ik charge assignment took %lf s.\n",t );

  t = MPI_Wtime();
  assign_forces( 1.0, s, &p, d_ik, f, 0);
  t -= MPI_Wtime();

  t_ik += t;

  printf("Ik force assignment took %lf s.\n",t );

  Free_data(d_ik);

  d_ad = Init_ad( s, &p );

  t = MPI_Wtime();
  assign_charge_and_derivatives( s, &p, d_ad, 0);
  t -= MPI_Wtime();

  t_ad += t;

  printf("Ad charge assignment took %lf s.\n",t );

  t = MPI_Wtime();
  assign_forces_ad( 1.0, s, &p, d_ad, f, 0);
  t -= MPI_Wtime();

  t_ad += t;

  printf("Ad force assignment took %lf s.\n",t );

  printf("Total time:\n");
  printf("Ik: %lf\n", -t_ik);
  printf("Ad: %lf\n", -t_ad);

  Free_data(d_ad);
}
