#include <math.h>
#include <string.h>

#include "common.h"

#include "realpart.h"


inline FLOAT_TYPE AS_erfc_part(double d)
{
#define AS_a1  0.254829592
#define AS_a2 -0.284496736
#define AS_a3  1.421413741
#define AS_a4 -1.453152027
#define AS_a5  1.061405429
#define AS_p   0.3275911
  double t;
  
  t = 1.0 / (1.0 + AS_p * d);
  
  return t * (AS_a1 + t * (AS_a2 + t * (AS_a3 + t * (AS_a4 + t * AS_a5) ) ) );
}


void Realteil( system_t *s, parameters_t *p, forces_t *f )
{
    /* Zwei Teilchennummern: */
    int t1,t2;
    /* Minimum-Image-Abstand: */
    FLOAT_TYPE dx,dy,dz,r;
    /* Staerke der elegktrostatischen Kraefte */
    FLOAT_TYPE fak;
    /* Zur Approximation der Fehlerfunktion */
    FLOAT_TYPE ar,erfc_teil;
    FLOAT_TYPE lengthi = 1.0/s->length;
    const FLOAT_TYPE wupi = 1.77245385090551602729816748334;

#pragma omp parallel for private(dx,dy,dz, r, ar, erfc_teil, fak, t2)
    for (t1=0; t1<s->nparticles; t1++) {
        for (t2=0; t2<s->nparticles; t2++) {
	  if(t1 == t2)
	    continue;

	  dx = s->p->x[t1] - s->p->x[t2];
	  dx -= round(dx*lengthi)*s->length;
	  dy = s->p->y[t1] - s->p->y[t2];
	  dy -= round(dy*lengthi)*s->length;
	  dz = s->p->z[t1] - s->p->z[t2];
	  dz -= round(dz*lengthi)*s->length;
	  
	  r = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
	  if (r<=p->rcut)
            {
	      ar= p->alpha*r;
	      //t = 1.0 / (1.0 + p*ar);
	      //erfc_teil = t*(a1+t*(a2+t*(a3+t*(a4+t*a5))));
	      erfc_teil = erfc(ar);
	      fak = s->q[t1]*s->q[t2]*
		(erfc_teil/r+(2.0*p->alpha/wupi)*exp(-ar*ar))/SQR(r);
	      
	      f->f_r->x[t1] += fak*dx;
	      f->f_r->y[t1] += fak*dy;
	      f->f_r->z[t1] += fak*dz;
            }
        }
    }
}

static inline void build_neighbor_list_for_particle(system_t *s, parameters_t *p, data_t *d, vector_array_t *buffer, int *neighbor_id_buffer, FLOAT_TYPE *charges_buffer, int id) {
    int i, j, np=0;
    FLOAT_TYPE r, dx, dy, dz;
    FLOAT_TYPE lengthi = 1.0/s->length;
    neighbor_list_t *neighbor_list = d->neighbor_list;

    for (i=id+1;i<s->nparticles;i++) {
        dx = s->p->x[id] - s->p->x[i];
        dx -= round(dx*lengthi)*s->length;
        dy = s->p->y[id] - s->p->y[i];
        dy -= round(dy*lengthi)*s->length;
        dz = s->p->z[id] - s->p->z[i];
        dz -= round(dz*lengthi)*s->length;

        r = sqrt(SQR(dx) + SQR(dy) + SQR(dz));

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
    int i;

    // Define and allocate buffers (nessecary due to unknow number of neigbors per particles).

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
 
    // Find neighbors for each particle.

    for (i=0;i<s->nparticles;i++)
      build_neighbor_list_for_particle( s, p, d, position_buffer, neighbor_id_buffer, charges_buffer, i );

    // NULL terminate the list

    neighbor_list[s->nparticles].n = -1;

    // Free buffers

    Free_vector_array(position_buffer);
    fftw_free(neighbor_id_buffer);
    fftw_free(charges_buffer);
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
    const FLOAT_TYPE wupii = 1.0/wupi;

    neighbor_list_t *neighbor_list = d->neighbor_list;

    //#pragma omp  parallel for private(dx, dy, dz, ar, r, erfc_teil, fak, j)

    for (i=0; i<s->nparticles-1; i++)
        for (j=0; j<neighbor_list[i].n; j++)
        {
            dx = s->p->fields[0][i] - neighbor_list[i].p->x[j];
            dx -= round(dx*lengthi)*s->length;
            dy = s->p->fields[1][i] - neighbor_list[i].p->y[j];
            dy -= round(dy*lengthi)*s->length;
            dz = s->p->fields[2][i] - neighbor_list[i].p->z[j];
            dz -= round(dz*lengthi)*s->length;

            r2 = SQR(dx) + SQR(dy) + SQR(dz);
            if (r2<=rcut2)
            {
	      r = sqrt( r2 );
	      ar= p->alpha*r;

	      //erfc_teil = erfc(ar);
              //  fak = s->q[i]*neighbor_list[i].q[j]*
	      //	  (erfc_teil/r+(2*p->alpha/wupi)*exp(-ar*ar))/ r2;

	      erfc_teil = AS_erfc_part( ar ) / r;
	      fak = s->q[i] * neighbor_list[i].q[j] * exp( - ar * ar ) * (erfc_teil + 2.0 * p->alpha * wupii ) / r2;

		//#pragma omp atomic
                f->f_r->x[i] += fak*dx;
		//#pragma omp atomic
                f->f_r->y[i] += fak*dy;
		//#pragma omp atomic
                f->f_r->z[i] += fak*dz;
		//#pragma omp atomic
                f->f_r->x[neighbor_list[i].id[j]] -= fak*dx;
		//#pragma omp atomic
                f->f_r->y[neighbor_list[i].id[j]] -= fak*dy;
		//#pragma omp atomic
                f->f_r->z[neighbor_list[i].id[j]] -= fak*dz;
            }
        }
}

void Free_neighborlist(data_t *d) {
  int i = 0;
  while(d->neighbor_list[i].n >= 0) {
    Free_vector_array(d->neighbor_list[i].p);
    fftw_free(d->neighbor_list[i].id);
    fftw_free(d->neighbor_list[i].q);
    i++;
  }
  fftw_free(d->neighbor_list);
}

FLOAT_TYPE Realspace_error( const system_t *s, const parameters_t *p )
{
  return (2.0*s->q2*exp(-SQR(p->rcut * p->alpha))) / (sqrt((double)s->nparticles* p->rcut ));
}
