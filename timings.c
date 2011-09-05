#include "timings.h"

#include <stdio.h>

#include <mpi/mpi.h>

#define MAX_TIMERS 10

struct {
  int size;
  double timers[MAX_TIMERS];
} timer_stack;

void start_timer(void) {
  if(timer_stack.size >= MAX_TIMERS) {
    fprintf(stderr, "Maximum number of timers (%d) reached!", MAX_TIMERS);
    return;
  }
  timer_stack.timers[timer_stack.size++] = MPI_Wtime();
}

double stop_timer(void) {
  if(timer_stack.size == 0)
    return -1.0;
  return MPI_Wtime() - timer_stack.timers[--timer_stack.size];
}

