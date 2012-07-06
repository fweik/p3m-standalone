#include <math.h>
#include <string.h>

#include "common.h"

#include "realpart.h"


int to_left=0;
int to_right=0;

static int *count_neighbors( system_t *s, parameters_t *p )
{
    /* Zwei Teilchennummern: */
    int t1,t2;
    /* Minimum-Image-Abstand: */
    FLOAT_TYPE dx,dy,dz,r;
    FLOAT_TYPE lengthi = 1.0/s->length;
    int *nb = fftw_malloc(s->nparticles*sizeof(int));

    memset(nb, 0, s->nparticles * sizeof(int));

    for (t1=0; t1<s->nparticles; t1++) {
        for (t2=0; t2<s->nparticles; t2++) {
	  if(t1 == t2)
	    continue;

	  dx = s->p->x[t1] - s->p->x[t2];
	  dx -= ROUND(dx*lengthi)*s->length;
	  dy = s->p->y[t1] - s->p->y[t2];
	  dy -= ROUND(dy*lengthi)*s->length;
	  dz = s->p->z[t1] - s->p->z[t2];
	  dz -= ROUND(dz*lengthi)*s->length;
	  
	  r = SQRT(SQR(dx) + SQR(dy) + SQR(dz));
	  if (r<=p->rcut)
            {
	      nb[t1]++;
            }
        }
    }
    return nb;
}


void Realteil( system_t *s, parameters_t *p, forces_t *f )
{
    /* Zwei Teilchennummern: */
    int t1,t2;
    /* Minimum-Image-Abstand: */
    FLOAT_TYPE dx,dy,dz,r, r2;
    /* Staerke der elegktrostatischen Kraefte */
    FLOAT_TYPE fak;
    /* Zur Approximation der Fehlerfunktion */
    FLOAT_TYPE ar,erfc_teil;
    FLOAT_TYPE lengthi = 1.0/s->length;
    const FLOAT_TYPE wupi = 1.77245385090551602729816748334;

    for (t1=0; t1<s->nparticles; t1++) {
        for (t2=0; t2<s->nparticles; t2++) {
	  if(t1 == t2)
	    continue;

	  dx = s->p->x[t1] - s->p->x[t2];
	  dx -= ROUND(dx*lengthi)*s->length;
	  dy = s->p->y[t1] - s->p->y[t2];
	  dy -= ROUND(dy*lengthi)*s->length;
	  dz = s->p->z[t1] - s->p->z[t2];
	  dz -= ROUND(dz*lengthi)*s->length;
	  
	  r = SQRT(SQR(dx) + SQR(dy) + SQR(dz));
	  if (r<=p->rcut)
            {
	      ar= p->alpha*r;

	      erfc_teil = ERFC(ar);
	      r2 = SQR(r);
	      fak = s->q[t1]*s->q[t2]*
		(erfc_teil/r+(2.0*p->alpha/wupi)*EXP(-ar*ar))/r2;
	      
	      f->f_r->x[t1] += fak*dx;
	      f->f_r->y[t1] += fak*dy;
	      f->f_r->z[t1] += fak*dz;

	      s->energy += 0.5 * s->q[t1] * s->q[t2] * erfc_teil / r;
            }
        }
    }
}

static void build_neighbor_list_for_particle(system_t *s, parameters_t *p, data_t *d, vector_array_t *buffer, int *neighbor_id_buffer, FLOAT_TYPE *charges_buffer, int id) {
  int i, j, np=0, last=id;
    FLOAT_TYPE r, dx, dy, dz;
    FLOAT_TYPE lengthi = 1.0/s->length;
    neighbor_list_t *neighbor_list = d->neighbor_list;

      to_left++;
      if(i < 0) {
	i += s->nparticles;
	if( i == id)
	  break;
      }
      dx = s->p->x[id] - s->p->x[i];
      dx -= ROUND(dx*lengthi)*s->length;
      if(dx > p->rcut) {
	break;
      }
      dy = s->p->y[id] - s->p->y[i];
      dy -= ROUND(dy*lengthi)*s->length;
      dz = s->p->z[id] - s->p->z[i];
      dz -= ROUND(dz*lengthi)*s->length;

      r = SQRT(SQR(dx) + SQR(dy) + SQR(dz));

      if (r<=p->rcut) {
	neighbor_id_buffer[np] = i;
	for (j=0;j<3;j++) {
	  buffer->fields[j][np] = s->p->fields[j][i];
	}
	charges_buffer[np] = s->q[i];
	last = i;
	np++;
      }
    }

    for (i=id+1;i!=id;i++) {
      to_right++;
      if(i >= s->nparticles) {
	i -= s->nparticles;
	if( i == id )
	  break;
      }
      if( i == last)
	break;
      dx = s->p->x[id] - s->p->x[i];
      dx -= ROUND(dx*lengthi)*s->length;
      if(fabs(dx) > p->rcut) {
	break;
      }
      dy = s->p->y[id] - s->p->y[i];
      dy -= ROUND(dy*lengthi)*s->length;
      dz = s->p->z[id] - s->p->z[i];
      dz -= ROUND(dz*lengthi)*s->length;

      r = SQRT(SQR(dx) + SQR(dy) + SQR(dz));

      if (r<=p->rcut) {
	neighbor_id_buffer[np] = i;
	for (j=0;j<3;j++) {
	  buffer->fields[j][np] = s->p->fields[j][i];
	}
	charges_buffer[np] = s->q[i];
	np++;
      }
    }
    neighbor_list[id].p = Init_vector_array(np);
    neighbor_list[id].q = Init_array(np, sizeof(FLOAT_TYPE));
    neighbor_list[id].id = Init_array(np, sizeof(int));

    for (j=0;j<3;j++)
        memcpy(neighbor_list[id].p->fields[j], buffer->fields[j], np*sizeof(FLOAT_TYPE));

    memcpy(neighbor_list[id].q, charges_buffer, np*sizeof(FLOAT_TYPE));
    memcpy(neighbor_list[id].id, neighbor_id_buffer, np*sizeof(int));

    neighbor_list[id].n = np;
}

void Init_neighborlist(system_t *s, parameters_t *p, data_t *d) {
  int i, *nb;

    // Define and allocate buffers (nessecary due to unknow number of neigbors per particle).

    int *neighbor_id_buffer = NULL;
    vector_array_t *position_buffer;
    FLOAT_TYPE *charges_buffer = NULL;
    neighbor_list_t *neighbor_list;

    neighbor_id_buffer = Init_array(s->nparticles, sizeof(int));
    position_buffer = Init_vector_array(s->nparticles);
    charges_buffer = Init_array(s->nparticles, sizeof(FLOAT_TYPE));

    // Allocate the actual list.

    d->neighbor_list = Init_array(s->nparticles + 1, sizeof(neighbor_list_t));
    neighbor_list = d->neighbor_list;   

    // Sort particles

    sort_particles(s);

    // Find neighbors for each particle.

    to_left = to_right = 0;

    for (i=0;i<s->nparticles;i++)
      build_neighbor_list_for_particle( s, p, d, position_buffer, neighbor_id_buffer, charges_buffer, i );

    // NULL terminate the list

    neighbor_list[s->nparticles].n = -1;

    //    printf("to_left %d to_right %d\n", to_left, to_right);

    // Free buffers

    Free_vector_array(position_buffer);
    FFTW_FREE(neighbor_id_buffer);
    FFTW_FREE(charges_buffer);
}

void Realpart_neighborlist(system_t *s, parameters_t *p, data_t *d, forces_t *f )
{
    int i,j;
    /* Minimum-Image-Abstand: */
    FLOAT_TYPE dx,dy,dz,r, r2, rcut2 = SQR ( p->rcut );
    /* Staerke der elegktrostatischen Kraefte */
    FLOAT_TYPE fak;
    /* Zur Approximation der Fehlerfunktion */
    FLOAT_TYPE ar,erfc_teil;
    FLOAT_TYPE lengthi = 1.0/s->length;
    const FLOAT_TYPE wupi = 1.77245385090551602729816748334;

    neighbor_list_t *neighbor_list = d->neighbor_list;

    for (i=0; i<s->nparticles; i++)
        for (j=0; j<neighbor_list[i].n; j++)
        {
            dx = s->p->fields[0][i] - neighbor_list[i].p->x[j];
            dx -= ROUND(dx*lengthi)*s->length;
            dy = s->p->fields[1][i] - neighbor_list[i].p->y[j];
            dy -= ROUND(dy*lengthi)*s->length;
            dz = s->p->fields[2][i] - neighbor_list[i].p->z[j];
            dz -= ROUND(dz*lengthi)*s->length;

            r = SQRT(SQR(dx) + SQR(dy) + SQR(dz));
            if ( r <= p->rcut )
            {
	      ar= p->alpha*r;

	      erfc_teil = ERFC(ar);
	      r2 = SQR(r);
	      fak = s->q[i]*neighbor_list[i].q[j]*
		(erfc_teil/r+(2.0*p->alpha/wupi)*EXP(-ar*ar))/r2;
	      
	      f->f_r->x[i] += fak*dx;
	      f->f_r->y[i] += fak*dy;
	      f->f_r->z[i] += fak*dz;
            }
        }
}

void Free_neighborlist(data_t *d) {
  int i = 0;
  while(d->neighbor_list[i].n >= 0) {
    Free_vector_array(d->neighbor_list[i].p);
    FFTW_FREE(d->neighbor_list[i].id);
    FFTW_FREE(d->neighbor_list[i].q);
    i++;
  }
  FFTW_FREE(d->neighbor_list);
}

FLOAT_TYPE Realspace_error( const system_t *s, const parameters_t *p )
{
  return (2.0*s->q2*EXP(-SQR(p->rcut * p->alpha))) / (SQRT((double)s->nparticles* p->rcut*s->length*s->length*s->length ));
}
