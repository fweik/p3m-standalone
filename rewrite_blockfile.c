#include <stdio.h>
#include <string.h>

#define FLOAT_TYPE double

#define MAXL 512

int main(int argc, char **argv) {
  FILE *b, *p, *f;
  char buf[MAXL];
  int id;
  FLOAT_TYPE fx, fy, fz;
  FLOAT_TYPE x,y,z;
  FLOAT_TYPE Q;

  if(argc != 4) {
    fprintf(stderr, "usage: %s <blockfile> <positionsfile> <forcesfile>\n", argv[0]);
    return 128;
  }

  b = fopen( argv[1], "r");
  p = fopen( argv[2], "w");
  f = fopen( argv[3], "w");

  do {
    fgets(buf, MAXL, b);
      } while(strncmp("{part", buf, 5) != 0);

  while(!feof(b)) {
    int Tl = 0;
    fgets(buf, MAXL, b);
    Tl = sscanf(buf, " { %d %lf %lf %lf %lf %lf %lf %lf } ", &id, &x, &y, &z, &Q, &fx, &fy, &fz);
    if(Tl == 8) {
      fprintf(p, "%.22f %.22f %.22f %.22f\n", x,y,z,Q);
      fprintf(f, "10.0 %.22f %.22f %.22f\n", fx, fy, fz);
    }
  }
  fclose(b);
  fclose(p);
  fclose(f);
}
