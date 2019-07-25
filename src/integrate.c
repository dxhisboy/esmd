#include <data.h>
#include <multiproc.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <pair_lj.h>
#include <thermo.h>
#include <loops.h>
#include <timer.h>
//#define DEBUG_THIS_FILE
#include <log.h>
void initial_integrate_nve(esmd_t *md){
  timer_start("init_integrate");
  box_t *box = md->box;
  //integrate_conf_t *integrate_conf = &(md->integrate_conf);
  potential_conf_t *pot_conf = md->pot_conf;
  areal dt = md->dt;
  areal dtf = dt * 0.5;
  areal *rmass = pot_conf->rmass;
  double totalv = 0;
  for (int kk = 0; kk < box->nlocal[2]; kk ++){
    for (int jj = 0; jj < box->nlocal[1]; jj ++){
      for (int ii = 0; ii < box->nlocal[0]; ii ++){

        int cell_off = get_cell_off(box, ii, jj, kk);
        cell_t *cell = box->cells + cell_off;
        celldata_t *data = box->celldata + cell_off;
        areal (*v)[3] = data->v;
        areal (*x)[3] = data->x;
        areal (*f)[3] = data->f;
        int *type = data->type;
        for (int i = 0; i < cell->natoms; i ++){
          totalv += fabs(v[i][0]) + fabs(v[i][1]) + fabs(v[i][2]);
          v[i][0] += dtf * f[i][0] * rmass[type[i]];
          v[i][1] += dtf * f[i][1] * rmass[type[i]];
          v[i][2] += dtf * f[i][2] * rmass[type[i]];
          x[i][0] += dt * v[i][0];
          x[i][1] += dt * v[i][1];
          x[i][2] += dt * v[i][2];
        }
      }
    }
  }
  timer_stop("init_integrate");
  //debug("totalv = %g\n", totalv);
}

void final_integrate_nve(esmd_t *md){
  timer_start("final_integrate");
  box_t *box = md->box;
  //integrate_conf_t *integrate_conf = &(md->integrate_conf);
  potential_conf_t *pot_conf = md->pot_conf;
  areal dt = md->dt;
  areal dtf = dt * 0.5;
  areal *rmass = pot_conf->rmass;
  double totalv = 0;
  for (int kk = 0; kk < box->nlocal[2]; kk ++){
    for (int jj = 0; jj < box->nlocal[1]; jj ++){
      for (int ii = 0; ii < box->nlocal[0]; ii ++){
	int cell_off = get_cell_off(box, ii, jj, kk);
	cell_t *cell = box->cells + cell_off;
	celldata_t *data = box->celldata + cell_off;
	areal (*v)[3] = data->v;
	areal (*f)[3] = data->f;
	int *type = data->type;
	for (int i = 0; i < cell->natoms; i ++){
	  v[i][0] += dtf * f[i][0] * rmass[type[i]];
	  v[i][1] += dtf * f[i][1] * rmass[type[i]];
	  v[i][2] += dtf * f[i][2] * rmass[type[i]];
	  totalv += fabs(v[i][0]) + fabs(v[i][1]) + fabs(v[i][2]);
	}
      }
    }
  }
  debug("totalv = %g\n", totalv);
  timer_stop("final_integrate");
}

#include <data.h>
void esmd_export_atoms(esmd_t *md){
  timer_start("export_atoms");
  box_t *box = md->box;
  areal *rlcell = box->rlcell;
  //integrate_conf_t *integrate_conf = &(md->integrate_conf);
  potential_conf_t *pot_conf = md->pot_conf;
  areal dt = md->dt;
  areal dtf = dt * 0.5;
  areal *rmass = pot_conf->rmass;
  int total_export = 0;
  for (int kk = 0; kk < box->nlocal[2]; kk ++){
    for (int jj = 0; jj < box->nlocal[1]; jj ++){
      for (int ii = 0; ii < box->nlocal[0]; ii ++){

	int cell_off = get_cell_off(box, ii, jj, kk);
	cell_t *cell = box->cells + cell_off;
	celldata_t *data = box->celldata + cell_off;
	areal (*v)[3] = data->v;
	areal (*x)[3] = data->x;
	areal (*f)[3] = data->f;
	areal *q = data->q;
	int *export = data->export;
	int *type = data->type;
	areal bbox[2][3];
	bbox[0][0] = cell->bbox_ideal[0][0] - TINY;
	bbox[0][1] = cell->bbox_ideal[0][1] + TINY;
	bbox[0][2] = cell->bbox_ideal[0][2] - TINY;
	bbox[1][0] = cell->bbox_ideal[1][0] + TINY;
	bbox[1][1] = cell->bbox_ideal[1][1] - TINY;
	bbox[1][2] = cell->bbox_ideal[1][2] + TINY;
	
	int natoms_new = 0, export_ptr = CELL_SIZE;
	for (int i = 0; i < cell->natoms; i ++){
	  int ecode = 0;
	  if (x[i][0] < bbox[0][0]) ecode -= 1;
	  else if (x[i][0] > bbox[1][0]) ecode += 1;
	  if (x[i][1] < bbox[0][1]) ecode -= 3;
	  else if (x[i][1] > bbox[1][1]) ecode += 3;
	  if (x[i][2] < bbox[0][2]) ecode -= 9;
	  else if (x[i][2] > bbox[1][2]) ecode += 9;

	  int i_new;
	  if (ecode != 0) {
	    export_ptr --;
	    i_new = export_ptr;
	  } else {
	    i_new = natoms_new;
	    natoms_new ++;
	  }
	  x[i_new][0] = x[i][0];
	  x[i_new][1] = x[i][1];
	  x[i_new][2] = x[i][2];
	  f[i_new][0] = f[i][0];
	  f[i_new][1] = f[i][1];
	  f[i_new][2] = f[i][2];
	  v[i_new][0] = v[i][0];
	  v[i_new][1] = v[i][1];
	  v[i_new][2] = v[i][2];
	  type[i_new] = type[i];
	  q[i_new] = q[i];
	  export[i_new] = ecode;
	}

	total_export += CELL_SIZE - export_ptr;
	cell->natoms = natoms_new;
	cell->export_ptr = export_ptr;
      }
    }
  }
  timer_stop("export_atoms");
  debug("%d\n", total_export);
}

inline int esmd_import_atoms_code(box_t *box, int self_off, int neigh_off, int icode) {
  cell_t *cell_self = box->cells + self_off;
  celldata_t *data_self = box->celldata + self_off;
  cell_t *cell_neigh = box->cells + neigh_off;
  celldata_t *data_neigh = box->celldata + neigh_off;
  areal *rlcell = box->rlcell;
  areal (*v)[3] = data_self->v;
  areal (*x)[3] = data_self->x;
  areal (*f)[3] = data_self->f;
  areal *q = data_self->q;
  int *type = data_self->type;

  areal (*ve)[3] = data_neigh->v;
  areal (*xe)[3] = data_neigh->x;
  areal (*fe)[3] = data_neigh->f;
  areal *qe = data_neigh->q;
  int *typee = data_neigh->type;
  int *export = data_neigh->export;
  int natoms_new = cell_self->natoms;
  for (int i = cell_neigh->export_ptr; i < CELL_SIZE; i ++){
    if (export[i] == icode){
      x[natoms_new][0] = xe[i][0];
      x[natoms_new][1] = xe[i][1];
      x[natoms_new][2] = xe[i][2];
      f[natoms_new][0] = fe[i][0];
      f[natoms_new][1] = fe[i][1];
      f[natoms_new][2] = fe[i][2];
      v[natoms_new][0] = ve[i][0];
      v[natoms_new][1] = ve[i][1];
      v[natoms_new][2] = ve[i][2];
      type[natoms_new] = typee[i];
      q[natoms_new] = qe[i];
      natoms_new ++;
    }
  }
  int nimport = natoms_new - cell_self->natoms;
  cell_self->natoms = natoms_new;
  return nimport;
}

void esmd_import_atoms(esmd_t *md){
  timer_start("import_atoms");
  box_t *box = md->box;
  ESMD_CELL_ITER(box, {
      for (int dz = -1; dz <= 1; dz ++) {
	for (int dy = -1; dy <= 1; dy ++){
	  for (int dx = -1; dx <= 1; dx ++){
	    int icode = dx * (-1) + dy * (-3) + dz * (-9);
	    int neighoff = get_cell_off(box, ii + dx, jj + dy, kk + dz);
            if ((1 << icode + 13) & box->cells[neighoff].emask)
              esmd_import_atoms_code(box, celloff, neighoff, icode);
	  }
	}
      }
    });
  timer_stop("import_atoms");
}

void integrate(esmd_t *md) {
  if (md->step % md->nthermo == 0){
    timer_start("compute");
    compute_kinetic_local(md);
    esmd_global_accumulate(md);
    thermo_compute(md);
    timer_stop("compute");
    master_info("step: %d eng: %f temp: %f press: %f\n", md->step, md->thermo.eng, md->thermo.temp, md->thermo.press);
  }
  md->accu_local.virial = 0;
  md->accu_local.epot = 0;
  md->accu_local.kinetic = 0;
  
  initial_integrate_nve_sw(md);
  
  esmd_export_atoms_sw(md);
  
  esmd_exchange_cell(md, LOCAL_TO_HALO, CELL_META | CELL_X | CELL_T | CELL_V | CELL_E, TRANS_ADJ_X | TRANS_EXPORTS);
  
  esmd_import_atoms(md);
  
  esmd_exchange_cell(md, LOCAL_TO_HALO, CELL_META | CELL_X | CELL_T | CELL_V, TRANS_ADJ_X | TRANS_ATOMS);
  int evflag = 0;
  if ((md->step + 1) % md->nthermo == 0) evflag = 3;
  pair_lj_force_sw(md, evflag);

  esmd_exchange_cell(md, HALO_TO_LOCAL, CELL_F, TRANS_INC_F | TRANS_ATOMS);
  
  //esmd_exchange_cell(md, LOCAL_TO_HALO, CELL_F);
  final_integrate_nve_sw(md);
  md->step ++;
}
