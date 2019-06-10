#include <data.h>
#include <math.h>
#include <stdlib.h>

void esmd_cell_setup_global(esmd_t *esmd, areal x, areal y, areal z){
  box_t *box = &(esmd->box);
  box->lglobal[0] = x;
  box->lglobal[1] = y;
  box->lglobal[2] = z;

  ireal lcell = box->lcell = esmd->pair_conf.cutoff / (NCELL_CUT + NCELL_SKIN);
  ireal rlcell = box->rlcell = 1 / box->lcell;
  
  box->nglobal[0] = (int)ceil(box->lglobal[0] * rlcell - TINY);
  box->nglobal[1] = (int)ceil(box->lglobal[1] * rlcell - TINY);
  box->nglobal[2] = (int)ceil(box->lglobal[2] * rlcell - TINY);
}

void esmd_cell_setup_local(esmd_t *esmd){

  box_t *box = &(esmd->box);
  
  int halo = NCELL_CUT;

  int nallx = box->nall[0], nally = box->nall[1], nallz = box->nall[2];
  int ncells = box->nall[0] * box->nall[1] * box->nall[2];

  
  cell_data_t *cell_data = (cell_data_t*)malloc(sizeof(cell_data_t) * ncells);
  
  cell_t *cells = (cell_t*)malloc(sizeof(cell_data_t) * ncells);
  
  esmd->box.cells = cells;

  areal lcell = box->lcell, rlcell = box->rlcell;

  for (int ii = 0; ii < nallx; ii ++) {
    for (int jj = 0; jj < nally; jj ++){
      for (int kk = 0; kk < nallz; kk ++){
        int i = ii - halo, j = jj - halo, k = kk - halo;
        cell_t *cell = esmd->box.cells + (ii * nally + jj) * nallz + kk;
        cell->natoms = 0;
        cell->nreplicas = 0;
        cell->bbox_ideal[0][0] = i * lcell;
        cell->bbox_ideal[0][1] = j * lcell;
        cell->bbox_ideal[0][2] = k * lcell;
        cell->bbox_ideal[1][0] = (i + 1) * lcell;
        cell->bbox_ideal[1][1] = (j + 1) * lcell;
        cell->bbox_ideal[1][2] = (k + 1) * lcell;
        cell->data = cell_data + (ii * nally + jj) * nallz + kk;
      }
    }
  }
}
