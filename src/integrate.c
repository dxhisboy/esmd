#include <data.h>
#include <multiproc.h>
#include <math.h>
#include <assert.h>
void initial_integrate_nve(esmd_t *md){
  box_t *box = &(md->box);
  integrate_conf_t *integrate_conf = &(md->integrate_conf);
  pair_conf_t *pair_conf = &(md->pair_conf);
  areal dt = integrate_conf->dt;
  areal dtf = dt * 0.5;
  areal *rmass = pair_conf->rmass;
  double totalv = 0;

  for (int ii = 0; ii < box->nlocal[0]; ii ++)
    for (int jj = 0; jj < box->nlocal[1]; jj ++)
      for (int kk = 0; kk < box->nlocal[2]; kk ++){
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
    printf("totalv = %g\n", totalv);
}

void final_integrate_nve(esmd_t *md){
  box_t *box = &(md->box);
  integrate_conf_t *integrate_conf = &(md->integrate_conf);
  pair_conf_t *pair_conf = &(md->pair_conf);
  areal dt = integrate_conf->dt;
  areal dtf = dt * 0.5;
  areal *rmass = pair_conf->rmass;
  double totalv = 0;
  for (int ii = 0; ii < box->nlocal[0]; ii ++)
    for (int jj = 0; jj < box->nlocal[1]; jj ++)
      for (int kk = 0; kk < box->nlocal[2]; kk ++){
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
  printf("totalv = %g\n", totalv);
}

#include <data.h>
void esmd_export_atoms(esmd_t *md){
  box_t *box = &(md->box);
  areal *rlcell = box->rlcell;
  integrate_conf_t *integrate_conf = &(md->integrate_conf);
  pair_conf_t *pair_conf = &(md->pair_conf);
  areal dt = integrate_conf->dt;
  areal dtf = dt * 0.5;
  areal *rmass = pair_conf->rmass;
  int total_export = 0;
  for (int ii = 0; ii < box->nlocal[0]; ii ++)
    for (int jj = 0; jj < box->nlocal[1]; jj ++)
      for (int kk = 0; kk < box->nlocal[2]; kk ++){
        int cell_off = get_cell_off(box, ii, jj, kk);
        cell_t *cell = box->cells + cell_off;
        celldata_t *data = box->celldata + cell_off;
        areal (*v)[3] = data->v;
        areal (*x)[3] = data->x;
        areal (*f)[3] = data->f;
        areal *q = data->q;
        int *export = data->export;
        int *type = data->type;

        int natoms_new = 0, export_ptr = CELL_SIZE;
        for (int i = 0; i < cell->natoms; i ++){
          /* if (x[i][0] < 0) x[i][0] += box->lglobal[0]; */
          /* if (x[i][1] < 0) x[i][1] += box->lglobal[1]; */
          /* if (x[i][2] < 0) x[i][2] += box->lglobal[2]; */
          /* if (x[i][0] > box->lglobal[0]) x[i][0] -= box->lglobal[0]; */
          /* if (x[i][1] > box->lglobal[1]) x[i][1] -= box->lglobal[1]; */
          /* if (x[i][2] > box->lglobal[2]) x[i][2] -= box->lglobal[2]; */
          int cellx = floor(x[i][0] * rlcell[0]);
          int celly = floor(x[i][1] * rlcell[1]);
          int cellz = floor(x[i][2] * rlcell[2]);
          if (cellx >= box->nglobal[0]) cellx -= box->nglobal[0];
          if (celly >= box->nglobal[1]) celly -= box->nglobal[1];
          if (cellz >= box->nglobal[2]) cellz -= box->nglobal[2];
          if (cellx < 0) cellx += box->nglobal[0];
          if (celly < 0) celly += box->nglobal[1];
          if (cellz < 0) cellz += box->nglobal[2];
          int new_cell_off = get_cell_off(box, cellx, celly, cellz);
          int i_new;
          if (cell_off != new_cell_off){
            export_ptr --;
            i_new = export_ptr;
            export[i_new] = new_cell_off;
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
        }
        total_export += CELL_SIZE - export_ptr;
        cell->natoms = natoms_new;
        cell->export_ptr = export_ptr;
        assert(export_ptr >= natoms_new);
      }
  printf("%d\n", total_export);
}

void esmd_import_atoms(esmd_t *md){
  box_t *box = &(md->box);
  areal *rlcell = box->rlcell;
  integrate_conf_t *integrate_conf = &(md->integrate_conf);
  pair_conf_t *pair_conf = &(md->pair_conf);
  areal dt = integrate_conf->dt;
  areal dtf = dt * 0.5;
  areal *rmass = pair_conf->rmass;
  int total_import = 0;
  for (int ii = 0; ii < box->nlocal[0]; ii ++)
    for (int jj = 0; jj < box->nlocal[1]; jj ++)
      for (int kk = 0; kk < box->nlocal[2]; kk ++){
        int cell_off = get_cell_off(box, ii, jj, kk);
        cell_t *cell = box->cells + cell_off;
        celldata_t *data = box->celldata + cell_off;
        areal (*v)[3] = data->v;
        areal (*x)[3] = data->x;
        areal (*f)[3] = data->f;
        areal *q = data->q;
        int *type = data->type;

        int natoms_new = cell->natoms;
        for (int dx = -1; dx <= 1; dx ++)
          for (int dy = -1; dy <= 1; dy ++)
            for (int dz = -1; dz <= 1; dz ++) {
              int neigh_off = get_cell_off(box, ii + dx, jj + dy, kk + dz);
              cell_t *cell_neigh = box->cells + neigh_off;
              celldata_t *data_neigh = box->celldata + neigh_off;
              areal (*ve)[3] = data_neigh->v;
              areal (*xe)[3] = data_neigh->x;
              areal (*fe)[3] = data_neigh->f;
              areal *qe = data_neigh->q;
              int *typee = data_neigh->type;
              int *export = data_neigh->export;
              for (int i = cell_neigh->export_ptr; i < CELL_SIZE; i ++){
                if (export[i] == cell_off){
                  x[natoms_new][0] = xe[i][0];
                  x[natoms_new][1] = xe[i][1];
                  x[natoms_new][2] = xe[i][2];
                  int safe = (x[natoms_new][0] + TINY >= cell->bbox_ideal[0][0] && x[natoms_new][0] -TINY <= cell->bbox_ideal[1][0] &&
                         x[natoms_new][1] + TINY >= cell->bbox_ideal[0][1] && x[natoms_new][1] -TINY <= cell->bbox_ideal[1][1] &&
                         x[natoms_new][2] + TINY >= cell->bbox_ideal[0][2] && x[natoms_new][2] -TINY <= cell->bbox_ideal[1][2]);
                  if (!safe){
                    printf("%d %d %d unsafely imported %f %f %f\n", ii, jj, kk, x[natoms_new][0], x[natoms_new][1], x[natoms_new][2]);
                  }
                  f[natoms_new][0] = fe[i][0];
                  f[natoms_new][1] = fe[i][1];
                  f[natoms_new][2] = fe[i][2];
                  v[natoms_new][0] = ve[i][0];
                  v[natoms_new][1] = ve[i][1];
                  v[natoms_new][2] = ve[i][2];
                  type[natoms_new] = typee[i];
                  q[natoms_new] = qe[i];
                  natoms_new ++;
                  total_import ++;
                }
              }
            }
        cell->natoms = natoms_new;
      }
  printf("%d\n", total_import);
}

void integrate(esmd_t *md) {
  initial_integrate_nve(md);
  esmd_export_atoms(md);
  esmd_exchange_cell(md, LOCAL_TO_HALO, CELL_META | CELL_X | CELL_T | CELL_E | CELL_V, TRANS_ADJ_X);
  esmd_import_atoms(md);
  esmd_exchange_cell(md, LOCAL_TO_HALO, CELL_META | CELL_X | CELL_T | CELL_E | CELL_V, TRANS_ADJ_X);
  pair_lj_force(md);
  esmd_exchange_cell(md, HALO_TO_LOCAL, CELL_F, TRANS_INC_F);
  final_integrate_nve(md);
}
