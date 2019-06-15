#include <data.h>
#include <math.h>
#include <stdlib.h>

void esmd_box_setup_global(esmd_t *md, areal x, areal y, areal z){
  box_t *box = &(md->box);
  box->lglobal[0] = x;
  box->lglobal[1] = y;
  box->lglobal[2] = z;

  ireal lcell = box->lcell = md->pair_conf.cutoff / (NCELL_CUT + NCELL_SKIN);
  ireal rlcell = box->rlcell = 1 / box->lcell;
  
  box->nglobal[0] = (int)ceil(box->lglobal[0] * rlcell + TINY);
  box->nglobal[1] = (int)ceil(box->lglobal[1] * rlcell + TINY);
  box->nglobal[2] = (int)ceil(box->lglobal[2] * rlcell + TINY);
}


inline void min(areal a, areal b) {
  if (a < b) return a;
  return b;
}
void esmd_box_setup_local(esmd_t *md){

  box_t *box = &(md->box);
  
  int halo = NCELL_CUT;

  int nallx = box->nall[0], nally = box->nall[1], nallz = box->nall[2];
  int ncells = box->nall[0] * box->nall[1] * box->nall[2];
  int nlocalx = box->nlocal[0], nlocaly = box->nlocal[1], nlocalz = box->nlocal[2];
  
  celldata_t *celldata = (celldata_t*)malloc(sizeof(celldata_t) * ncells);
  
  cell_t *cells = (cell_t*)malloc(sizeof(celldata_t) * ncells);
  
  box->cells = cells;
  box->celldata = celldata;
  areal lcell = box->lcell, rlcell = box->rlcell;
  areal lx = box->lglobal[0], ly = box->lglobal[1], lz = box->lglobal[2];
  for (int ii = 0; ii < nlocalx; ii ++) {
    for (int jj = 0; jj < nlocaly; jj ++) {
      for (int kk = 0; kk < nlocalz; kk ++){
        int cell_off = get_cell_off(box, ii, jj, kk);
        cell_t *cell = box->cells + cell_off;
        cell->natoms = 0;
        cell->nreplicas = 0;
        cell->bbox_ideal[0][0] = i * lcell;
        cell->bbox_ideal[0][1] = j * lcell;
        cell->bbox_ideal[0][2] = k * lcell;
        cell->bbox_ideal[1][0] = (i + 1) * lcell;
        cell->bbox_ideal[1][1] = (j + 1) * lcell;
        cell->bbox_ideal[1][2] = (k + 1) * lcell;        
      }
    }
  }
  /* for (int ii = 0; ii < nallx; ii ++) { */
  /*   for (int jj = 0; jj < nally; jj ++){ */
  /*     for (int kk = 0; kk < nallz; kk ++){ */
  /*       int i = ii - halo, j = jj - halo, k = kk - halo; */
  /*       cell_t *cell = md->box.cells + (ii * nally + jj) * nallz + kk; */
  /*       cell->natoms = 0; */
  /*       cell->nreplicas = 0; */
  /*       cell->bbox_ideal[0][0] = i * lcell; */
  /*       cell->bbox_ideal[0][1] = j * lcell; */
  /*       cell->bbox_ideal[0][2] = k * lcell; */
  /*       cell->bbox_ideal[1][0] = (i + 1) * lcell; */
  /*       cell->bbox_ideal[1][1] = (j + 1) * lcell; */
  /*       cell->bbox_ideal[1][2] = (k + 1) * lcell; */
  /*       cell->data = celldata + (ii * nally + jj) * nallz + kk; */
  /*     } */
  /*   } */
  /* } */
  
}
