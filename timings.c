/**    Copyright (C) 2011,2012,2013 Florian Weik <fweik@icp.uni-stuttgart.de>

       This program is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

       This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this program.  If not, see <http://www.gnu.org/licenses/>. **/

#include "timings.h"

#include <stdio.h>

#include <mpi.h>

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

