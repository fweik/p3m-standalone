
#include <time.h>

double wtime(void) {
  double wtime;
  struct timeval tv;
  gettimeofday(&tv, NULL);
  wtime = tv.tv_sec;
  wtime += (double)tv.tv_usec / 1000000.0;
  return wtime;
}
