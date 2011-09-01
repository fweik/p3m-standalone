#include <math.h>
#include <string.h>

#include "realpart.h"

static neighbor_list_t *neighbor_list;

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
    for (t1=0; t1<s->nparticles-1; t1++) {
        for (t2=t1+1; t2<s->nparticles; t2++) {
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

#pragma omp atomic
                f->f_r->x[t1] += fak*dx;
#pragma omp atomic
                f->f_r->y[t1] += fak*dy;
#pragma omp atomic
                f->f_r->z[t1] += fak*dz;
#pragma omp atomic
                f->f_r->x[t2] -= fak*dx;
#pragma omp atomic
                f->f_r->y[t2] -= fak*dy;
#pragma omp atomic
                f->f_r->z[t2] -= fak*dz;

            }
        }
    }
}

static inline void build_neighbor_list_for_particle(system_t *s, parameters_t *p, vector_array_t *buffer, int *neighbor_id_buffer, FLOAT_TYPE *charges_buffer, int id) {
    int i, j, np=0;
    FLOAT_TYPE r, dx, dy, dz;
    FLOAT_TYPE lengthi = 1.0/s->length;

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

void Init_neighborlist(system_t *s, parameters_t *p) {
    int i;

    // Define and allocate buffers (nessecary due to unknow number of neigbors per particles).

    int *neighbor_id_buffer = NULL;
    vector_array_t *position_buffer;
    FLOAT_TYPE *charges_buffer = NULL;

    neighbor_id_buffer = Init_array(s->nparticles, sizeof(int));
    position_buffer = Init_vector_array(s->nparticles);
    charges_buffer = Init_array(s->nparticles, sizeof(FLOAT_TYPE));

    // Allocate the actual list.

    neighbor_list = Init_array(s->nparticles, sizeof(neighbor_list_t));

    // Find neighbors for each particle.

    for (i=0;i<s->nparticles;i++)
        build_neighbor_list_for_particle( s, p, position_buffer, neighbor_id_buffer, charges_buffer, i );

    // Free buffers

    Free_vector_array(position_buffer);
    free(neighbor_id_buffer);
    free(charges_buffer);
}

void Realpart_neighborlist(system_t *s, parameters_t *p, forces_t *f )
{
    int i,j;
    /* Minimum-Image-Abstand: */
    FLOAT_TYPE dx,dy,dz,r;
    /* Staerke der elegktrostatischen Kraefte */
    FLOAT_TYPE fak;
    /* Zur Approximation der Fehlerfunktion */
    FLOAT_TYPE ar,erfc_teil;
    FLOAT_TYPE lengthi = 1.0/s->length;
    const FLOAT_TYPE wupi = 1.77245385090551602729816748334;

#pragma omp parallel for private(dx, dy, dz, ar, r, erfc_teil, fak, j)
    FLOAT_TYPE t;
    for (i=0; i<s->nparticles-1; i++)
        for (j=0; j<neighbor_list[i].n; j++)
        {
            dx = s->p->x[i] - neighbor_list[i].p->x[j];
            dx -= round(dx*lengthi)*s->length;
            dy = s->p->y[i] - neighbor_list[i].p->y[j];
            dy -= round(dy*lengthi)*s->length;
            dz = s->p->z[i] - neighbor_list[i].p->z[j];
            dz -= round(dz*lengthi)*s->length;

            r = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
            if (r<=p->rcut)
            {
                ar= p->alpha*r;
                //t = 1.0 / (1.0 + p*ar);
                erfc_teil = t*(a1+t*(a2+t*(a3+t*(a4+t*a5))));
                //erfc_teil = erfc(ar);
                fak = s->q[i]*neighbor_list[i].q[j]*
                      (erfc_teil/r+(2*p->alpha/wupi)*exp(-ar*ar))/SQR(r);

#pragma omp atomic
                f->f_r->x[i] += fak*dx;
#pragma omp atomic
                f->f_r->y[i] += fak*dy;
#pragma omp atomic
                f->f_r->z[i] += fak*dz;
#pragma omp atomic
                f->f_r->x[neighbor_list[i].id[j]] -= fak*dx;
#pragma omp atomic
                f->f_r->y[neighbor_list[i].id[j]] -= fak*dy;
#pragma omp atomic
                f->f_r->z[neighbor_list[i].id[j]] -= fak*dz;
            }
        }
}
