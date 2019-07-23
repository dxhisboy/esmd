#include <math.h>
#include <data.h>
#include <loops.h>
#include <multiproc.h>
#include <log.h>
double temperature(esmd_t *md){
  box_t *box = md->box;
  areal mv2_tot = 0;
  ESMD_CELL_ITER(box, {
      for (int i = 0; i < cell->natoms; i ++){
	ireal mass = md->pot_conf->mass[type[i]];
	mv2_tot += v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2] * mass;
      }
    });
  areal mv2_gbl;
  esmd_global_sum_scalar(md, &mv2_gbl, mv2_tot);
  return mv2_gbl / (md->natoms * 3 - 3);
}

double compute_kinetic_local(esmd_t *md){
  box_t *box = md->box;
  areal mv2_tot = 0;
  ESMD_CELL_ITER(box, {
      for (int i = 0; i < cell->natoms; i ++){
	ireal mass = md->pot_conf->mass[type[i]];
	mv2_tot += v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2] * mass;
      }
    });
  md->accu_local.kinetic = mv2_tot * 0.5;
}

void thermo_init(esmd_t *md){
  areal *lglobal = md->box->lglobal;
  if (md->utype == UNIT_LJ){
    md->thermo.t_scale = 1.0 / (md->natoms * 3 - 3);
    md->thermo.p_scale = 1.0 / 3 / (lglobal[0] * lglobal[1] * lglobal[2]);
    md->thermo.e_scale = 0.5;
    md->thermo.dof_boltz = (md->natoms * 3 - 3);
  }
}
void thermo_compute(esmd_t *md){
  thermo_t *thermo = &(md->thermo);
  accumulate_t *accu = &(md->accu_global);
  thermo->temp = accu->kinetic * 2 * thermo->t_scale;
  thermo->press = (thermo->temp * thermo->dof_boltz + accu->virial) * thermo->p_scale;
  thermo->eng = accu->epot / md->natoms * thermo->e_scale;
}

void scale_to_temp(esmd_t *md, areal t_req){
  box_t *box = md->box;
  areal t0 = temperature(md);
  ireal scale = sqrt(t_req / t0);
  ESMD_CELL_ITER(box, {
      for (int i = 0; i < cell->natoms; i ++){
	v[i][0] *= scale;
	v[i][1] *= scale;
	v[i][2] *= scale;
	if (x[i][0] == 0 && x[i][1] == 0 && x[i][2] == 0){
	  debug("%f %f %f\n", v[i][0], v[i][1], v[i][2]);
	}	
      }
    });
}
