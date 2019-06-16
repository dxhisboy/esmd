#ifndef DATA_H_
#define DATA_H_

#include <cppdefs.h>
#include <mpi.h>
typedef AREAL areal;
typedef IREAL ireal;
enum cell_fields {
  CELL_META = 1,
  CELL_X = 2,
  CELL_V = 4,
  CELL_Q = 8,
  CELL_T = 16
};

typedef struct celldata {
  areal f[CELL_SIZE][3];
  areal v[CELL_SIZE][3];
  areal x[CELL_SIZE][3];
  int type[CELL_SIZE];
#ifdef CHARGED
  ireal q[CELL_SIZE];
#endif
} celldata_t;

#define CELL_DATA_XTQ_SIZE (sizeof(celldata_t) - ((celldata_t *)NULL).x)

typedef struct cellforce {
  areal f[CELL_SIZE][3];
} cellforce_t;

typedef struct cell {
  areal bbox_ideal[2][3], bbox_real[2][3];
  areal trans[3];
  int natoms;
  int nreplicas;
  cellforce_t **replicas;
  celldata_t *data;
} cell_t;

typedef struct box {
  int nlocal[3], nall[3], nglobal[3];
  int offset[3];
  areal lcell[3], rlcell[3];
  areal lglobal[3]; //, llocal[3], lall[3];
  //areal olocal[3], oall[3];
  cell_t *cells;
  celldata_t *celldata;
} box_t;
typedef ireal type_table[MAX_TYPES];
typedef ireal pair_table[MAX_TYPES][MAX_TYPES];
typedef ireal angle_table[MAX_TYPES][MAX_TYPES][MAX_TYPES];
typedef ireal dihedral_table[MAX_TYPES][MAX_TYPES][MAX_TYPES][MAX_TYPES];

typedef struct pair_conf {
  ireal cutoff;
  struct {
    ireal coef6, coef12;
  } lj_param;
} pair_conf_t;

typedef struct multiproc {
  int npx, npy, npz, np;
  int pidx, pidy, pidz, pid;
  MPI_Comm comm;
} multiproc_t;

#include <memory.h>
typedef struct esmd {
  mempool_t force_pool;
  pair_conf_t pair_conf;
  box_t box;
  multiproc_t mpp;
} esmd_t;


#define get_cell_off(boxptr, i, j, k) ((((i) + NCELL_CUT) * (boxptr)->nall[1] \
                                        + (j) + NCELL_CUT) * (boxptr)->nall[2] \
                                       + (k) + NCELL_CUT)
#endif
