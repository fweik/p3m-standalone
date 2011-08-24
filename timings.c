#include "timings.h"

#include <mpi/mpi.h>

double start_time;

void start_timer(void) {
  start_time = MPI_Wtime();  
}

double stop_timer(void) {
  return MPI_Wtime() - start_time; 
}

