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
