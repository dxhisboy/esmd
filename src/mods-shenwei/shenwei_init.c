#ifdef MPE
#include <athread.h>
#include <stdio.h>
#include <mpi.h>
#include <data.h>
#include <swdata.h>
#define LWPF_UNITS U(PAIR_LJ)
#undef inline
#include <lwpf2/lwpf2.h>
#define DEBUG_THIS_FILE
#include <log.h>
static inline void part1d(int n, int np, int pid, int *start, int *count){
  int pncell = n / np;
  int rncell = n % np;
  if (pid < rncell) {
    *start = (pncell + 1) * pid;
    *count = pncell + 1;
  } else {
    *start = pncell * pid + rncell;
    *count = pncell;
  }
}

perf_config_t conf;
void shenwei_init(esmd_t *md){
  athread_init();
  conf.pcrc = PCRC_ALL;
  conf.pcr0 = PC0_CNT_INST;
  conf.pcr1 = PC1_CYCLE;
  conf.pcr2 = PC2_CNT_GLD;
  lwpf_init(&conf);
  box_t *box = md->box;
  int *nlocal = box->nlocal;
  int *nall = box->nall;
  int ncellsall = nall[0] * nall[1] * nall[2];
  swdata_t *swdata = esmd_malloc(sizeof(swdata_t), "shenwei data");
  //partition the box
  for (int k = 0; k < 4; k ++){
    for (int j = 0; j < 4; j ++){
      for (int i = 0; i < 4; i ++){
        int ipe = k * 16 + j * 4 + i;
        part1d(nlocal[2], 4, k, swdata->start[ipe] + 2, swdata->count[ipe] + 2);
        part1d(nlocal[1], 4, j, swdata->start[ipe] + 1, swdata->count[ipe] + 1);
        part1d(nlocal[0], 4, i, swdata->start[ipe] + 0, swdata->count[ipe] + 0);
      }
    }
  }
  
  for (int i = 0; i < ncellsall; i ++){
    box->cells[i].pemask = 0;
  }
  int nupdate = 0;
  for (int ipe = 0; ipe < NCPE; ipe ++){
    int *start = swdata->start[ipe];
    int *count = swdata->count[ipe];
    int ks = start[2] - NCELL_CUT;
    int js = start[1] - NCELL_CUT;
    int is = start[0] - NCELL_CUT;
    int ke = start[2] + count[2] +  NCELL_CUT;
    int je = start[1] + count[1] +  NCELL_CUT;
    int ie = start[0] + count[0] +  NCELL_CUT;
    unsigned long long ipemask = 1ULL << ipe;
    for (int kk = ks; kk < ke; kk ++){
      for (int jj = js; jj < je; jj ++){
        for (int ii = is; ii < ie; ii ++){
          int celloff = get_cell_off(box, ii, jj, kk);
          cell_t *cell = box->cells + celloff;
          if ((cell->pemask & ipemask) == 0) nupdate ++;
          //swdata->pemask[celloff] |= 1ULL << ipe;
          cell->pemask |= ipemask;
        }
      }
    }
  }
  //debug("%d %d %d %d %d\n", ncellsall, nupdate, swdata->count[0][0], swdata->count[0][1], swdata->count[0][2]);
  //mempool_init(&(swdata->fpool), sizeof(areal) * CELL_SIZE * 3, nupdate, "force pool");
  swdata->frep = esmd_malloc(sizeof(areal) * CELL_SIZE * 3 * nupdate, "force replicas");
  int rep_head = 0;
  for (int i = 0; i < ncellsall; i ++){
    cell_t *cell = box->cells + i;
    int nrep;
    asm("ctpop %1, %0\n\t" : "=r"(nrep) : "r"(cell->pemask));
    cell->nreplicas = nrep;
    cell->frep = swdata->frep + rep_head;
    rep_head ++;
  }
  swdata->nupdates = nupdates;
}

void shenwei_destroy(){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0){
    lwpf_report_summary(stdout, &conf);
  }
}
#endif
