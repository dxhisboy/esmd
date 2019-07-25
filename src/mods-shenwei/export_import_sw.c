#include <data.h>
#ifdef CPE
#include <slave.h>
#include <dma_macros.h>
#include <cal.h>
void esmd_export_atoms_cpe(esmd_t *gl_md){
  /* esmd_t md; */
  /* pe_get(gl_md, &md, sizeof(esmd_t)); */
  /* dma_syn(); */
  dma_init();
  box_t box;
  pe_get(gl_md->box, &box, sizeof(box_t));
  dma_syn();
  areal *rlcell = box.rlcell;

  areal x[CELL_SIZE][3], v[CELL_SIZE][3], f[CELL_SIZE][3];
  ireal q[CELL_SIZE];
  int type[CELL_SIZE], export[CELL_SIZE];
  int *nlocal = box.nlocal;
  int ncells = nlocal[0] * nlocal[1] * nlocal[2];
  for (int icell = _MYID; icell < ncells; icell += 64){
    int kk = icell / (nlocal[0] * nlocal[1]);
    int jj = icell / nlocal[0] % nlocal[1];
    int ii = icell % nlocal[0];

    /* for (int kk = _ROW; kk < box.nlocal[2]; kk += 8){ */
    /*   for (int jj = _COL; jj < box.nlocal[1]; jj += 8){ */
    /*     for (int ii = 0; ii < box.nlocal[0]; ii ++){ */
    int cell_off = get_cell_off((&box), ii, jj, kk);
    /* cell_t *cell = box.cells + cell_off; */
    celldata_t *data = box.celldata + cell_off;
    /* areal (*v)[3] = data->v; */
    /* areal (*x)[3] = data->x; */
    /* areal (*f)[3] = data->f; */
    /* areal *q = data->q; */
    /* int *export = data->export; */
    /* int *type = data->type; */
    cell_t cell;
    pe_get(box.cells + cell_off, &cell, sizeof(cell_t));
    dma_syn();
    pe_get(data->x, x, cell.natoms * sizeof(areal) * 3);
    pe_get(data->v, v, cell.natoms * sizeof(areal) * 3);
    pe_get(data->q, q, cell.natoms * sizeof(ireal));
    pe_get(data->type, type, cell.natoms * sizeof(int));
    dma_syn();
    areal bbox[2][3];
    bbox[0][0] = cell.bbox_ideal[0][0] - TINY;
    bbox[0][1] = cell.bbox_ideal[0][1] + TINY;
    bbox[0][2] = cell.bbox_ideal[0][2] - TINY;
    bbox[1][0] = cell.bbox_ideal[1][0] + TINY;
    bbox[1][1] = cell.bbox_ideal[1][1] - TINY;
    bbox[1][2] = cell.bbox_ideal[1][2] + TINY;
    int emask = 0;
    int natoms_new = 0, export_ptr = CELL_SIZE;
    for (int i = 0; i < cell.natoms; i ++){
      int ecode = 0;
      if (x[i][0] < bbox[0][0]) ecode -= 1;
      else if (x[i][0] > bbox[1][0]) ecode += 1;
      if (x[i][1] < bbox[0][1]) ecode -= 3;
      else if (x[i][1] > bbox[1][1]) ecode += 3;
      if (x[i][2] < bbox[0][2]) ecode -= 9;
      else if (x[i][2] > bbox[1][2]) ecode += 9;
      emask |= 1 << (ecode + 13);
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
      v[i_new][0] = v[i][0];
      v[i_new][1] = v[i][1];
      v[i_new][2] = v[i][2];
      type[i_new] = type[i];
      q[i_new] = q[i];
      export[i_new] = ecode;
    }
    int nexport = CELL_SIZE - export_ptr;
    cell.natoms = natoms_new;
    cell.export_ptr = export_ptr;
    cell.emask = emask;
    pe_put(box.cells + cell_off, &cell, sizeof(cell_t));

    if (nexport > 0){
      pe_put(data->x, x, sizeof(areal) * 3 * CELL_SIZE);
      pe_put(data->v, v, sizeof(areal) * 3 * CELL_SIZE);
      pe_put(data->q, q, sizeof(ireal) * CELL_SIZE);
      pe_put(data->type, type, sizeof(int) * CELL_SIZE);
      pe_put(data->export, export, sizeof(int) * CELL_SIZE);
    }
    dma_syn();
  }
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

void esmd_import_atoms_cpe(esmd_t *md){
  box_t *box = md->box;
  int *nlocal = box->nlocal;
  int ncells = nlocal[0] * nlocal[1] * nlocal[2];
  for (int icell = _MYID; icell < ncells; icell += 64){
    int kk = icell / (nlocal[0] * nlocal[1]);
    int jj = icell / nlocal[0] % nlocal[1];
    int ii = icell % nlocal[0];
    //ESMD_CELL_ITER(box, {
    int celloff = get_cell_off(box, ii, jj, kk);
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
  }
}

#endif

#ifdef MPE
#include <athread.h>
#include <timer.h>
extern void slave_esmd_export_atoms_cpe(esmd_t *md);
void esmd_export_atoms_sw(esmd_t *md){
  timer_start("export_atoms_sw");
  athread_spawn(esmd_export_atoms_cpe, md);
  athread_join();
  timer_stop("export_atoms_sw");
}
#endif
