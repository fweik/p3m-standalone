#include "tuning.h"

#include "timings.h"
#include "interpol.h"

parameters_t *Tune( method_t *m, system_t *s, FLOAT_TYPE precision ) {
  parameters_t *p = Init_array( 1, sizeof( parameters_t ) );
  // Parameter iteraters, best parameter set
  parameters_t it, p_best;
  forces_t *f = Init_forces(s->nparticles);

  data_t *d = NULL;

  FLOAT_TYPE rcut_max = 0.5 * s->length;
  FLOAT_TYPE rcut_step = RCUT_STEP * s->length;
  FLOAT_TYPE last_success = -1.0;

  FLOAT_TYPE best_time=1e250, time;

  int success = 0;
  int direction=1;

  for(it.mesh = MESH_MIN; it.mesh <= MESH_MAX; it.mesh++ ) {
    for(it.cao = CAO_MAX; it.cao >= CAO_MIN; it.cao--) {
      it.cao3 = it.cao * it.cao * it.cao;
      it.ip = it.cao - 1;

      Free_data(d);
      d = m->Init( s, &it );

      last_success = -1.0;

      for(it.rcut = 0; rcut_step >= RCUT_STEP_MIN * s->length; it.rcut += direction * rcut_step ) {
	if(it.rcut > rcut_max)
	  break;
        if(last_success >= 0.0)
	  rcut_step /= 2.0;
	if( m->Error( s, &it ) <= precision ) {
	  last_success = it.rcut;
	  direction = -1;
	}
	else
	  direction = 1;
      }
      if(last_success >= 0.0) {
	it.rcut = last_success;   
	success = 1;

	start_timer();
	m->Kspace_force( s, p, d, f );
	time = stop_timer();

	if( time < best_time ) {
	  p_best.mesh = it.mesh;
	  p_best.rcut = it.rcut;
	  p_best.cao = it.cao;
	}
      } else {
	// If we didn't get the precision here, there is no point going to lower cao.
	break;
      }
    }
    if( time > 1.2*best_time ) {
      // No further gain expected with bigger mesh;
      break;
    }
  }

  Free_data(d);
  Free_forces(f);

  if( success == 0 )
    return NULL;

  *p = p_best;
  p->ip = p->cao - 1;
  p->cao3 = p->cao * p->cao * p->cao;

  return p;
}  


