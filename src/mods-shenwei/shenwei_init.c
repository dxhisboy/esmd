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
  int ncells = nlocal[0] * nlocal[1] * nlocal[2];
  swdata_t *swdata = esmd_malloc(sizeof(swdata_t), "shenwei data");
  md->platformdata = swdata;
  //partition the box
  /* for (int k = 0; k < 8; k ++){ */
  /*   for (int j = 0; j < 4; j ++){ */
  /*     for (int i = 0; i < 2; i ++){ */
  /*       int ipe = k * 8 + j * 2 + i; */
  /*       part1d(nlocal[2], 8, k, swdata->start[ipe] + 2, swdata->count[ipe] + 2); */
  /*       part1d(nlocal[1], 4, j, swdata->start[ipe] + 1, swdata->count[ipe] + 1); */
  /*       part1d(nlocal[0], 2, i, swdata->start[ipe] + 0, swdata->count[ipe] + 0); */
  /*     } */
  /*   } */
  /* } */
  
  for (int i = 0; i < ncellsall; i ++){
    box->cells[i].pemask = 0;
  }
  int nupdate = 0;

  for (int ipe = 0; ipe < NCPE; ipe ++){
    unsigned long long ipemask = 1ULL << ipe;
    for (int icellst = ipe * BLK; icellst < ncells; icellst += NCPE * BLK){
      int icelled = icellst + BLK;
      if (icelled > ncells) icelled = ncells;
      for (int icell = icellst; icell < icelled; icell ++){
        int kk = icell / (nlocal[0] * nlocal[1]);
        int jj = icell / nlocal[0] % nlocal[1];
        int ii = icell % nlocal[0];
	for (int dz = -NCELL_CUT; dz <= 0; dz ++) {
	  int dytop = (dz == 0) ? 0 : NCELL_CUT;
	  for (int dy = -NCELL_CUT; dy <= dytop; dy ++) {
	    int dxtop = (dz == 0 && dy == 0) ? 0 : NCELL_CUT;
	    for (int dx = -NCELL_CUT; dx <= dxtop; dx ++) {
              int celloff = get_cell_off(box, ii + dx, jj + dy, kk + dz);
              cell_t *cell = box->cells + celloff;
              if ((cell->pemask & ipemask) == 0) nupdate ++;
              //swdata->pemask[celloff] |= 1ULL << ipe;
              cell->pemask |= ipemask;
            }
          }
        }
      }
    }
  }
  //debug("%d %d %d %d %d\n", ncellsall, nupdate, swdata->count[0][0], swdata->count[0][1], swdata->count[0][2]);
  //mempool_init(&(swdata->fpool), sizeof(areal) * CELL_SIZE * 3, nupdate, "force pool");
  swdata->frep = esmd_malloc(sizeof(areal) * CELL_SIZE * 3 * nupdate, "force replicas");
  memset(swdata->frep, -1, sizeof(areal) * CELL_SIZE * 3 * nupdate);
  int rep_head = 0;
  for (int i = 0; i < ncellsall; i ++){
    cell_t *cell = box->cells + i;
    int nrep;
    asm("ctpop %1, %0\n\t" : "=r"(nrep) : "r"(cell->pemask));
    cell->nreplicas = nrep;
    cell->frep = swdata->frep + rep_head;
    rep_head += nrep;
  }
  //debug("%p %d\n", box->cells[5121].frep, box->cells[5121].nreplicas);
  swdata->nupdate = nupdate;
  int neigh0[NNEIGHBOR][3], neigh1[NNEIGHBOR][3];
  for (int dz = -NCELL_CUT; dz <= 0; dz ++) {
    int dytop = (dz == 0) ? 0 : NCELL_CUT;
    for (int dy = -NCELL_CUT; dy <= dytop; dy ++) {
      int dxtop = (dz == 0 && dy == 0) ? 0 : NCELL_CUT;
      for (int dx = -NCELL_CUT; dx <= dxtop; dx ++) {
	int off = dz * 9 + dy * 3 + dx + 13;
	neigh0[off][0] = 0 + dx;
	neigh0[off][1] = 0 + dy;
	neigh0[off][2] = 0 + dz;
	neigh1[off][0] = 1 + dx;
	neigh1[off][1] = 0 + dy;
	neigh1[off][2] = 0 + dz;
      }
    }
  }
  int *read_from = swdata->read_from, *store_to = swdata->store_to;
  int rwoff = 0;
  for (int i = 0; i < NNEIGHBOR; i ++) {
    read_from[i] = -1;
    store_to[i] = -1;
  }
  for (int i = 0; i < NNEIGHBOR; i ++) {
    for (int j = 0; j < NNEIGHBOR; j ++){
      if (neigh0[i][0] == neigh1[j][0] &&
	  neigh0[i][1] == neigh1[j][1] &&
	  neigh0[i][2] == neigh1[j][2]){
	store_to[i] = rwoff;
	read_from[j] = rwoff;
	debug("delta %d %d, rwoff %d\n", i, j, rwoff);
	rwoff ++;
      }
    }
  }
}

void shenwei_destroy(){
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0){
    lwpf_report_summary(stdout, &conf);
  }
}
#endif
