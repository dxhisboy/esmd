#include <math.h>
#include <data.h>
#include <loops.h>
#include <multiproc.h>
double temperature(esmd_t *md){
  box_t *box = &(md->box);
  areal mv2_tot = 0;
  ESMD_CELL_ITER(box, {
      for (int i = 0; i < cell->natoms; i ++){
	ireal mass = md->pair_conf.mass[type[i]];
	mv2_tot += v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2] * mass;
      }
    });
  areal mv2_gbl;
  esmd_global_sum_scalar(md, &mv2_gbl, mv2_tot);
  return mv2_gbl / (md->natoms * 3 - 3);
}

void scale_to_temp(esmd_t *md, areal t_req){
  box_t *box = &(md->box);
  areal t0 = temperature(md);
  ireal scale = sqrt(t_req / t0);
  ESMD_CELL_ITER(box, {
      for (int i = 0; i < cell->natoms; i ++){
	v[i][0] *= scale;
	v[i][1] *= scale;
	v[i][2] *= scale;
	if (x[i][0] == 0 && x[i][1] == 0 && x[i][2] == 0){
	  printf("%f %f %f\n", v[i][0], v[i][1], v[i][2]);
	}	
      }
    });
}
