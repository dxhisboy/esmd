#ifndef DATA_H_
#define DATA_H_

#include <camd_config.h>
typedef AREAL areal;
typedef IREAL ireal;

typedef struct cell_data {
  areal x[CELL_SIZE][3];
#ifdef CHARGED
  ireal q[CELL_SIZE];
#endif
#ifdef PER_ATOM_MASS
  ireal m[CELL_SIZE];
#endif
  int type[CELL_SIZE];
  int export[CELL_SIZE];
} cell_data_t;

typedef struct cell_force {
  areal f[CELL_SIZE][3];
} cell_force_t;

typedef struct cell {
  areal bbox_ideal[2][3], bbox_real[2][3];
  areal trans[3];
  int natoms;
  int nreplicas;
  cell_force_t **replicas;
  cell_data_t *data;
} cell_t;

typedef struct box {
  int nlocal[3], nall[3], nglobal[3];
  int offset[3];
  areal lcell, rlcell;
  areal lglobal[3], llocal[3], lall[3];
  areal olocal[3], oall[3];
  cell_t *cells;
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
#include <memory.h>

typedef struct camd {
  mempool_t force_pool;
  pair_conf_t pair_conf;
  box_t box;

} camd_t;


#endif
