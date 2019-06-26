#include <data.h>
#include <math.h>
#include <stdlib.h>

void esmd_box_setup_global(esmd_t *md, areal x, areal y, areal z){
  box_t *box = &(md->box);
  box->lglobal[0] = x;
  box->lglobal[1] = y;
  box->lglobal[2] = z;

  ireal lcell[3];
  
  ireal rlcell_max = (NCELL_CUT + NCELL_SKIN) / md->pair_conf.cutoff;
  
  box->nglobal[0] = (int)floor(box->lglobal[0] * rlcell_max + TINY);
  box->nglobal[1] = (int)floor(box->lglobal[1] * rlcell_max + TINY);
  box->nglobal[2] = (int)floor(box->lglobal[2] * rlcell_max + TINY);
  
  box->lcell[0] = box->lglobal[0] / box->nglobal[0];
  box->lcell[1] = box->lglobal[1] / box->nglobal[1];
  box->lcell[2] = box->lglobal[2] / box->nglobal[2];
  
  box->rlcell[0] = 1. / box->lcell[0];
  box->rlcell[1] = 1. / box->lcell[1];
  box->rlcell[2] = 1. / box->lcell[2];
  printf("lcell: %f %f %f\n", box->lcell[0], box->lcell[1], box->lcell[2]);
}


inline areal min(areal a, areal b) {
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
  areal *lcell = box->lcell, *rlcell = box->rlcell;

  areal lx = box->lglobal[0], ly = box->lglobal[1], lz = box->lglobal[2];
  for (int ii = 0; ii < nlocalx; ii ++) {
    for (int jj = 0; jj < nlocaly; jj ++) {
      for (int kk = 0; kk < nlocalz; kk ++){
        int cell_off = get_cell_off(box, ii, jj, kk);
        cell_t *cell = box->cells + cell_off;
        cell->natoms = 0;
        cell->nreplicas = 0;
        cell->bbox_ideal[0][0] = ii * lcell[0];
        cell->bbox_ideal[0][1] = jj * lcell[1];
        cell->bbox_ideal[0][2] = kk * lcell[2];
        cell->bbox_ideal[1][0] = (ii + 1) * lcell[0];
        cell->bbox_ideal[1][1] = (jj + 1) * lcell[1];
        cell->bbox_ideal[1][2] = (kk + 1) * lcell[2];
      }
    }
  }
}
