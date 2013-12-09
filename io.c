#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "types.h"
#include "io.h"
#include "common.h"

#include "tools/visit_writer.h"

void write_vtf(char *filename, system_t *s) {
  FILE *f = fopen( filename, "w");
  int i;

  if(f == NULL) {
    fprintf( stderr, "Could not open '%s' for writing.", filename);
    return;
  }

  for(i=0;i<s->nparticles;i++) {
    fprintf( f, "a %d t %c r 0.5\n", i, (s->q[i] < 0.0) ? 'O' : 'H');
  }

  fprintf( f, "timestep indexed\n");

  for(i=0;i<s->nparticles;i++) {
    fprintf( f, "%d %lf %lf %lf\n", i, FLOAT_CAST s->p->x[i], FLOAT_CAST s->p->y[i], FLOAT_CAST s->p->z[i] );
  }  
} 


void write_mesh(char *filename, FLOAT_TYPE *data, int *dims, FLOAT_TYPE *spacing, int data_size, const char *var_name) {
  int i,j,k,index_row,index_col;
  int centering[] = {1};
  float *data_col_maj=NULL;

  data_col_maj = malloc(dims[0]*dims[1]*dims[2]*data_size*sizeof(float));

  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++)
      for(i=0;i<dims[0];i++) {
	index_row = (dims[1]*dims[2]*i + dims[1]*j + k)*data_size;
	index_col = (dims[1]*dims[2]*k + dims[1]*j + i)*data_size;
	data_col_maj[index_col] =  data[index_row];
      }
  write_regular_mesh( filename, 0, dims, 1, &data_size, centering, &var_name, &data_col_maj);
}

void Read_exact_forces(system_t *s, char *filename)
{
    FILE *fp;
    int i, ret_val;
    double buf[6];
    int taw;

    printf("npart %d\n", s->nparticles);

    fp=fopen(filename, "r");

    if ((fp == NULL) || feof(fp)) {
        fprintf(stderr, "Could not open '%s' for reading.\n", filename);
        exit(127);
    }

    for (i=0; i<s->nparticles; i++) {
        ret_val = fscanf(fp,"%d\t%lf\t%lf\t%lf\n", &taw, buf, buf + 1, buf + 2);
	s->reference->f->x[i] = buf[0];
	s->reference->f->y[i] = buf[1];
	s->reference->f->z[i] = buf[2];
			 
	if((ret_val != 4)) {
          fprintf(stderr, "Error while reading file '%s' (%d) (%d)\n", filename, ret_val, i);
          exit(-1);
        }
    }
    fclose(fp);
}

system_t *Read_system(parameters_t *p, char *filename)
{
    /* Opens file 'filename' for reanding and reads system parameters,
     particle positions and charges. */

    FILE *fp;
    int i;

    double buf[4];
    FLOAT_TYPE Length;
    int n;

    int ret_val = 0;

    system_t *s;

    assert(p != NULL);
    assert(filename != NULL);

    fp=fopen(filename, "r");

    if ((fp == NULL) || feof(fp)) {
        fprintf(stderr, "Could not open '%s' for reading.\n", filename);
        exit(127);
    }

    //read system parameters
    ret_val += fscanf(fp,"# Teilchenzahl: %d\n",&n);
    ret_val += fscanf(fp,"# Len: %lf\n",buf);
    Length = buf[0];

    if(ret_val != 2) {
      fprintf(stderr, "Error while reading file '%s'\n", filename);
      exit(-1);
    }
    
    s = Init_system(n);
    s->length = Length;

    s->q2 = 0.0;
    /* Teilchenkoordinaten und -ladungen: */
    for (i=0; i<s->nparticles; i++) {
      ret_val = fscanf(fp,"%lf\t%lf\t%lf\t%lf\n", buf, buf + 1, buf + 2, buf + 3 );
      s->p->x[i] = buf[0];
      s->p->y[i] = buf[1];
      s->p->z[i] = buf[2];
      s->q[i] = buf[3];
	
        if(ret_val != 4) {
          fprintf(stderr, "Error while reading file '%s'\n", filename);
          exit(-1);
        }
        s->q2 += SQR(s->q[i]);
    }
    fclose(fp);

    return s;
}

void Write_exact_forces(system_t *s, char *forces_file) {
    FILE *fin;
    int i;

    fin = fopen(forces_file, "w");

    if (fin == NULL) {
        fprintf(stderr, "Could not open '%s' for writing!.\n", forces_file);
        exit(127);
    }

    for (i=0;i<s->nparticles;i++) {
        fprintf(fin, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                i, FLOAT_CAST s->reference->f->x[i], FLOAT_CAST s->reference->f->y[i], FLOAT_CAST s->reference->f->z[i],
		FLOAT_CAST s->reference->f_k->x[i], FLOAT_CAST s->reference->f_k->y[i], FLOAT_CAST s->reference->f_k->z[i]);
    }

    fclose(fin);
}

void Write_system_cuda( system_t *s, parameters_t *p, char *filename) {
    FILE *fp;

    assert(filename != NULL);

    fp=fopen(filename, "w");

    if ((fp == NULL) || feof(fp)) {
        fprintf(stderr, "Could not open '%s' for writing.\n", filename);
        exit(127);
    }
  
    fprintf(fp, "%d %d %d %lf %lf\n", s->nparticles, p->cao, p->mesh, p->alpha, s->length);

    for(int i = 0; i < s->nparticles; i++) {
      fprintf( fp, "%.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", FLOAT_CAST (s->p->x[i]), FLOAT_CAST (s->p->y[i]), FLOAT_CAST (s->p->z[i]), FLOAT_CAST (s->q[i]), FLOAT_CAST s->reference->f_k->x[i], FLOAT_CAST s->reference->f_k->y[i], FLOAT_CAST s->reference->f_k->z[i]);
    }
    fclose(fp);
}

void Write_system(system_t *s, char *filename)
{
    /* Opens file 'filename' for reanding and reads system parameters,
     particle positions and charges. */

    FILE *fp;
    int i;

    assert(filename != NULL);

    fp=fopen(filename, "w");

    if ((fp == NULL) || feof(fp)) {
        fprintf(stderr, "Could not open '%s' for writing.\n", filename);
        exit(127);
    }

    //read system parameters
    fprintf(fp,"# Teilchenzahl: %d\n", s->nparticles);
    fprintf(fp,"# Len: %lf\n", FLOAT_CAST  s->length);

    /* Teilchenkoordinaten und -ladungen: */
    for (i=0; i<s->nparticles; i++) {
        fprintf(fp,"%lf\t%lf\t%lf\t%lf\n", FLOAT_CAST (s->p->x[i]), FLOAT_CAST (s->p->y[i]), FLOAT_CAST (s->p->z[i]), FLOAT_CAST (s->q[i]));
    }
    fclose(fp);
}
