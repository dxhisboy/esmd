#ifndef DATA_H_
#define DATA_H_

#include <cppdefs.h>
#ifdef NOMPI
#else
#include <mpi.h>
#endif
typedef AREAL areal;
typedef IREAL ireal;
enum cell_fields {
  CELL_META = 1,
  CELL_X = 2,
  CELL_V = 4,
  CELL_Q = 8,
  CELL_T = 16,
  CELL_F = 32,
  CELL_E = 64
};

enum unit_type {
  UNIT_LJ,
  UNIT_REAL
};

typedef struct celldata {
  int export[CELL_SIZE];
  areal f[CELL_SIZE][3];
  areal v[CELL_SIZE][3];
  areal x[CELL_SIZE][3];
  int type[CELL_SIZE];
  ireal q[CELL_SIZE];
} celldata_t;

#define struct_off(type, field) ((size_t)(((type *)NULL).field))
#define CELL_DATA_XT_SIZE (struct_off(celldata_t, q) - struct_off(celldata_t, x))
#define CELL_DATA_XTQ_SIZE (sizeof(celldata_t) - struct_off(celldata_t, x))

typedef struct cellforce {
  areal f[CELL_SIZE][3];
} cellforce_t;

typedef struct cell {
  areal bbox_ideal[2][3];
  int natoms, nexport, export_ptr;
  int nreplicas, emask, padding;
  cellforce_t **replicas;
  celldata_t *data;
} cell_t;

enum cell_commtype {
  CT_INNER = 0,
  CT_OUTER = 1,
  CT_HALO
};
typedef struct box {
  int nlocal[3], nall[3], nglobal[3];
  int offset[3];
  areal lcell[3], rlcell[3];
  areal lglobal[3];
  areal llocal[3], olocal[3]; //, oall[3];
  cell_t *cells;
  celldata_t *celldata;
  int *celltype, *cellowner;
} box_t;
typedef ireal type_table[MAX_TYPES];
typedef ireal pair_table[MAX_TYPES][MAX_TYPES];
typedef ireal angle_table[MAX_TYPES][MAX_TYPES][MAX_TYPES];
typedef ireal dihedral_table[MAX_TYPES][MAX_TYPES][MAX_TYPES][MAX_TYPES];

typedef struct lj_param {
  pair_table c6, c12, cutoff2, ec6, ec12;
} lj_param_t;

typedef struct potential_conf {
  ireal cutoff;
  union {
    lj_param_t lj;
  } param;
  type_table rmass, mass;
} potential_conf_t;

enum integrate_type {
  FIX_NVE,
  FIX_NVT,
  FIX_NPT
};

typedef struct integrate_conf {
  int type;
  areal dt;
} integrate_conf_t;

typedef struct halo {
  int off[3][2], len[3];
  int ncells;
  areal translation[3];
  void *recv_buf, *send_buf;
  int neighbor;
  int send_tag, recv_tag;
  MPI_Request send_req, recv_req;
  MPI_Status send_stat, recv_stat;
} halo_t;

typedef struct multiproc {
  int npx, npy, npz, np;
  int pidx, pidy, pidz, pid;
  halo_t halo[MAX_COMM];  // axis * prev/next
  MPI_Comm comm;
} multiproc_t;

enum lattice_type {
  LAT_FCC
};

typedef struct lattice {
  areal (*offset)[3];
  areal lx, ly, lz;
  int natoms;
} lattice_t;

typedef struct lattice_config {
  enum lattice_type type;
  areal dens, scale;
  int nx, ny, nz;
  int *atom_types;
} lattice_conf_t;

typedef struct accumulate {
  areal epot, virial, kinetic;
} accumulate_t;

typedef struct thermo {
  areal eng, temp, press;
  areal e_scale, t_scale, p_scale, dof_boltz;
} thermo_t;

#include <memory.h>
typedef struct esmd {
  potential_conf_t *pot_conf;
  box_t *box;
  multiproc_t *mpp;
  int fix_type;
  areal dt;
  //integrate_conf_t integrate_conf;
  lattice_conf_t lat_conf;
  enum unit_type utype;
  accumulate_t accu_local, accu_global;
  thermo_t thermo;
  int natoms, step, nthermo;
} esmd_t;


/* #define get_cell_off(boxptr, i, j, k) ((((i) + NCELL_CUT) * (boxptr)->nall[1] \ */
/*                                         + (j) + NCELL_CUT) * (boxptr)->nall[2] \ */
/*                                        + (k) + NCELL_CUT) */

#define get_cell_off(boxptr, i, j, k) ((((k) + NCELL_CUT) * (boxptr)->nall[1] \
                                        + (j) + NCELL_CUT) * (boxptr)->nall[0] \
                                       + (i) + NCELL_CUT)
#endif
