#include <stdio.h>
#include <data.h>
#include <math.h>
#include <stdlib.h>
#include <util.h>
//#define DEBUG_THIS_FILE
#include <loops.h>
#include <log.h>

void esmd_set_box_size(esmd_t *md, areal x, areal y, areal z){
  box_t *box = &(md->box);
  box->lglobal[0] = x;
  box->lglobal[1] = y;
  box->lglobal[2] = z;  
}
void esmd_box_setup_global(esmd_t *md){
  box_t *box = &(md->box);
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
}

inline int get_cell_type_1d(int idx, int nlocal){
  if (idx < 0 || idx >= nlocal) return CT_HALO;
  if (idx < NCELL_CUT || idx >= nlocal - NCELL_CUT) return CT_OUTER;
  return CT_INNER;

}

void esmd_box_setup_local(esmd_t *md){

  box_t *box = &(md->box);
  
  int halo = NCELL_CUT;

  int nallx = box->nall[0], nally = box->nall[1], nallz = box->nall[2];
  int ncells = box->nall[0] * box->nall[1] * box->nall[2];
  int nlocalx = box->nlocal[0], nlocaly = box->nlocal[1], nlocalz = box->nlocal[2];

  box->llocal[0] = box->nlocal[0] * box->lcell[0];
  box->llocal[1] = box->nlocal[1] * box->lcell[1];
  box->llocal[2] = box->nlocal[2] * box->lcell[2];
  
  box->olocal[0] = box->offset[0] * box->lcell[0];
  box->olocal[1] = box->offset[1] * box->lcell[1];
  box->olocal[2] = box->offset[2] * box->lcell[2];
  
  celldata_t *celldata = (celldata_t*)esmd_malloc(sizeof(celldata_t) * ncells, "celldata");
  
  cell_t *cells = (cell_t*)esmd_malloc(sizeof(cell_t) * ncells, "cellmeta");
  int *celltype = (int*)esmd_malloc(sizeof(int) * ncells, "celltype");
  box->cells = cells;
  box->celldata = celldata;
  box->celltype = celltype;
  areal *lcell = box->lcell, *rlcell = box->rlcell;

  areal lx = box->lglobal[0], ly = box->lglobal[1], lz = box->lglobal[2];

  for (int kk = 0; kk < nlocalz; kk ++){
    for (int jj = 0; jj < nlocaly; jj ++) {
      for (int ii = 0; ii < nlocalx; ii ++) {
        int cell_off = get_cell_off(box, ii, jj, kk);
        cell_t *cell = box->cells + cell_off;
        cell->natoms = 0;
        cell->nreplicas = 0;
        cell->bbox_ideal[0][0] = (box->offset[0] + ii) * lcell[0];
        cell->bbox_ideal[0][1] = (box->offset[1] + jj) * lcell[1];
        cell->bbox_ideal[0][2] = (box->offset[2] + kk) * lcell[2];
        cell->bbox_ideal[1][0] = (box->offset[0] + ii + 1) * lcell[0];
        cell->bbox_ideal[1][1] = (box->offset[1] + jj + 1) * lcell[1];
        cell->bbox_ideal[1][2] = (box->offset[2] + kk + 1) * lcell[2];
      }
    }
  }
  for (int kk = -NCELL_CUT; kk < box->nlocal[2] + NCELL_CUT; kk ++){
    int typez = get_cell_type_1d(kk, box->nlocal[2]);
    for (int jj = -NCELL_CUT; jj < box->nlocal[1] + NCELL_CUT; jj ++){
      int typey = get_cell_type_1d(jj, box->nlocal[1]);
      for (int ii = -NCELL_CUT; ii < box->nlocal[0] + NCELL_CUT; ii ++){
	int typex = get_cell_type_1d(ii, box->nlocal[0]);
        int cell_off = get_cell_off(box, ii, jj, kk);
        celltype[cell_off] = max(typex, max(typey, typez));
      }
    }
  }

}

void box_add_atom(box_t *box, areal *x, areal *v, ireal q, int type){
  areal *rlcell = box->rlcell;
  int ci = floor(x[0] * rlcell[0] + TINY) - box->offset[0];
  int cj = floor(x[1] * rlcell[1] + TINY) - box->offset[1];
  int ck = floor(x[2] * rlcell[2] + TINY) - box->offset[2];
  int celloff = get_cell_off(box, ci, cj, ck);
  cell_t *cell = box->cells + celloff;
  celldata_t *celldata = box->celldata + celloff;
  int curatom = cell->natoms;
  celldata->x[curatom][0] = x[0];
  celldata->x[curatom][1] = x[1];
  celldata->x[curatom][2] = x[2];
  celldata->v[curatom][0] = v[0];
  celldata->v[curatom][1] = v[1];
  celldata->v[curatom][2] = v[2];
  //to be filled
  celldata->type[curatom] = 0;
  celldata->q[curatom] = q;
  celldata->f[curatom][0] = 0;
  celldata->f[curatom][1] = 0;
  celldata->f[curatom][2] = 0;
  cell->natoms ++;
}

void report_cell_info(esmd_t *md){
  box_t *box = &(md->box);
  int nmax = 0, nmin = 0x7fffffff, nsum = 0;
  ESMD_CELL_ITER(box, {
      if (cell->natoms > nmax) nmax = cell->natoms;
      if (cell->natoms < nmin) nmin = cell->natoms;
      nsum += cell->natoms;
      //if (cell->natoms == 32) debug("%d %d %d %d %d\n", ii, jj, kk, cell->natoms, md->mpp.pid);
    });
  
  int gbl_nmax, gbl_nmin, gbl_nsum;
  MPI_Reduce(&nmax, &gbl_nmax, 1, MPI_INT, MPI_MAX, 0, md->mpp.comm);
  MPI_Reduce(&nmin, &gbl_nmin, 1, MPI_INT, MPI_MIN, 0, md->mpp.comm);
  MPI_Reduce(&nsum, &gbl_nsum, 1, MPI_INT, MPI_SUM, 0, md->mpp.comm);
  double avg = gbl_nsum * 1.0 / (box->nglobal[0] * box->nglobal[1] * box->nglobal[2]);
  if (md->mpp.pid == 0){
    info("number of atoms per cell: avg=%f, min=%d, max=%d\n", avg, gbl_nmin, gbl_nmax);
  }
}
