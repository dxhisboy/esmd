#include <mpi.h>
#include <data.h>
inline void part1d(int n, int np, int pid, int *start, int *count){
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
void esmd_multiproc_part_cart(esmd_t *md, int npx, int npy, int npz, int pid){
  box_t *box = &(md->box);
  int pidx = pid / (npy * npz);
  int pidy = (pid - pidx * npy * npz) / npz;
  int pidz = pid - (pidx * npy + pidy) * npz;
  int pcellx = box->nglobal[0] / npx;
  int pcelly = box->nglobal[1] / npy;
  int pcellz = box->nglobal[2] / npz;
  md->mpp.pidx = pidx;
  md->mpp.pidy = pidy;
  md->mpp.pidz = pidz;
  md->mpp.pid = pid;
  md->mpp.npx  = npx ;
  md->mpp.npy  = npy ;
  md->mpp.npz  = npz ;
  md->mpp.np   = npx * npy * npz;
  int stx, cntx, sty, cnty, stz, cntz;
  part1d(box->nglobal[0], npx, pidx, &stx, &cntx);
  part1d(box->nglobal[1], npy, pidy, &sty, &cnty);
  part1d(box->nglobal[2], npz, pidz, &stz, &cntz);
  /* box->llocal[0] = cntx * box->lcell; */
  /* box->llocal[1] = cnty * box->lcell; */
  /* box->llocal[2] = cntz * box->lcell; */
  /* box->lall[0] = (cntx + halo * 2) * box->lcell; */
  /* box->lall[1] = (cnty + halo * 2) * box->lcell; */
  /* box->lall[2] = (cntz + halo * 2) * box->lcell; */
  
  box->nall[0] = cntx + NCELL_CUT * 2;
  box->nall[1] = cnty + NCELL_CUT * 2;
  box->nall[2] = cntz + NCELL_CUT * 2;
  
  box->nlocal[0] = cntx;
  box->nlocal[1] = cnty;
  box->nlocal[2] = cntz;

  box->offset[0] = stx;
  box->offset[1] = sty;
  box->offset[2] = stz;
}

inline int proc_3d_to_flat(multiproc_t *mpp, int pidx, int pidy, int pidz){
  return pidx * mpp->npy * mpp->npz + pidy * mpp->npz + pidz;
}

inline int check_offset();
#define TAG_ZN 30
#define TAG_ZP 31
#define TAG_YN 20
#define TAG_YP 21
#define TAG_XN 10
#define TAG_XP 11
void esmd_exchange_cell(esmd_t *md, int fields) {
  box_t *box = &(md->box);
  multiproc_t *mpp = &(md->mpp);
  int *nlocal = box->nlocal;
  int *nall = box->nall;
  int buffer_size = max(nall[0] * nall[1], nall[0] * nall[2], nall[1] * nall[2]);
  int entry_size = esmd_fields_size(fields);
  void *send_buf = esmd_malloc(buffer_size * entry_size * 2, "exchange_buffer");
  void *recv_buf = esmd_malloc(buffer_size * entry_size * 2, "exchange_buffer");
  MPI_Request send_req[2], recv_req[2];
  MPI_Request send_stat[2], recv_stat[2];
  size_t ncomm_z = entry_size * box->nloal[0] * box->nlocal[1];
  MPI_Irecv(recv_buf, ncomm_z, MPI_CHAR, proc_zn, TAG_ZN, mpp->comm, recv_req);
  MPI_Irecv(recv_buf + buffer_size, ncomm_z, MPI_CHAR, proc_zp, TAG_ZP, mpp->comm, recv_req + 1);
  
  for (int i = 0; i < box->nlocal[0]; i ++)
    for (int j = 0; j < box->nlocal[1]; j ++) {
      size_t bufoff_n = (i * nlocal[1] + j) * entry_size;
      size_t bufoff_p = bufoff_n + buffer_size;
      esmd_export_cell(md, send_buf + bufoff_n, fields, get_cell_off(box, i, j, 0));
      esmd_export_cell(md, send_buf + bufoff_p, fields, get_cell_off(box, i, j, nlocal[2] - 1));
    }
  int procz_zn, procz_zp, proc_zn, proc_zp;
  procz_zn = (mpp->pidz == 0) ? mpp->npz - 1 : mpp->pidz - 1;
  procz_zp = (mpp->pidz == mpp->npz - 1) ? 0 : mpp->pidz + 1;
  proc_zn = procz_zn + (mpp->pidx * mpp->npy + mpp->pidy) * mpp->npz;
  proc_zp = procz_zp + (mpp->pidx * mpp->npy + mpp->pidy) * mpp->npz;
  MPI_Isend(send_buf, ncomm_z, MPI_CHAR, proc_zn, TAG_ZP, mpp->comm, send_req);
  MPI_Isend(send_buf + buffer_size, ncomm_z, MPI_CHAR, proc_zn, TAG_ZP, mpp->comm, send_req + 1);

  for (int i = 0; i < box->nlocal[0]; i ++)
    for (int j = 0; j < box->nlocal[0]; j ++){
      size_t bufoff_n = (i * nlocal[1] + j) * entry_size;
      size_t bufoff_p = bufoff_n + buffer_size;
      esmd_import_cell(md, recv_buf + bufoff_n, fields, get_cell_off(box, i, j, -1));
      esmd_import_cell(md, recv_buf + bufoff_p, fields, get_cell_off(box, i, j, nlocal[2]));
    }
}
