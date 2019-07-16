#ifndef MULTIPROC_H_
#define MULTIPROC_H_
#include <data.h>
#include <util.h>
enum exchange_direction {
  LOCAL_TO_HALO = 0,
  HALO_TO_LOCAL = 1
};

//function signatures
void init_halo_unordered(halo_t *halo, esmd_t *md, int dx, int dy, int dz);
void init_comm_unordered(esmd_t *md);
void init_comm_ordered(esmd_t *md);
void esmd_multiproc_part_cart(esmd_t *md, int npx, int npy, int npz, int pid);
void esmd_comm_start(esmd_t *md, MPI_Comm comm, halo_t *halo, int dir, int fields, int flags);
void esmd_comm_finish(esmd_t *md, halo_t *halo, int dir, int fields, int flags);
void esmd_exchange_cell(esmd_t *md, int direction, int fields, int flags);
void esmd_exchange_cell_ordered(esmd_t *md, int direction, int fields, int flags);
void esmd_global_sum_vec(esmd_t *md, areal *result, areal *localvec);
void esmd_global_sum_scalar(esmd_t *md, areal *result, areal localval);
void esmd_global_accumulate(esmd_t *md);
//end function signatures
#endif
