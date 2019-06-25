#include <mpi.h>
#include <data.h>
#include <util.h>
#include <stdio.h>
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

#define TAG_ZN 30
#define TAG_ZP 31
#define TAG_YN 20
#define TAG_YP 21
#define TAG_XN 10
#define TAG_XP 11

static inline int get_neighbor(esmd_t *md, int dx, int dy, int dz, areal *off) {
  multiproc_t *mpp = &(md->mpp);
  int neigh_x = mpp->pidx + dx;
  int neigh_y = mpp->pidy + dy;
  int neigh_z = mpp->pidz + dz;
  off[0] = off[1] = off[2] = 0;
  if (neigh_x < 0) {
    neigh_x += mpp->npx;
    off[0] = -md->box.lglobal[0];
  }
  if (neigh_x >= mpp->npx) {
    neigh_x -= mpp->npx;
    off[0] = md->box.lglobal[0];
  }
  if (neigh_y < 0) {
    neigh_y += mpp->npy;
    off[1] = -md->box.lglobal[1];
  }
  if (neigh_y >= mpp->npy) {
    neigh_y -= mpp->npy;
    off[1] = md->box.lglobal[1];
  }
  if (neigh_z < 0) {
    neigh_z += mpp->npz;
    off[2] = -md->box.lglobal[2];
  }
  if (neigh_z >= mpp->npz) {
    neigh_z -= mpp->npz;
    off[2] = md->box.lglobal[2];
  }
  return (neigh_x * mpp->npy + neigh_y) * mpp->npz + neigh_z;
}
void esmd_exchange_cell(esmd_t *md, int fields) {
  box_t *box = &(md->box);
  multiproc_t *mpp = &(md->mpp);
  int *nlocal = box->nlocal;
  int *nall = box->nall;
  int buffer_size = max(nall[0] * nall[1], max(nall[0] * nall[2], nall[1] * nall[2])) * NCELL_CUT;
  size_t entry_size = esmd_fields_size(fields);
  void *send_buf_n = esmd_malloc(buffer_size * entry_size, "exchange_buffer");
  void *send_buf_p = esmd_malloc(buffer_size * entry_size, "exchange_buffer");
  void *recv_buf_n = esmd_malloc(buffer_size * entry_size, "exchange_buffer");
  void *recv_buf_p = esmd_malloc(buffer_size * entry_size, "exchange_buffer");
  
  MPI_Request send_req_n, send_req_p, recv_req_n, recv_req_p;
  MPI_Status send_stat_n, send_stat_p, recv_stat_n, recv_stat_p;

  areal off_n[3], off_p[3];
  
  size_t ncomm_z = entry_size * nlocal[0] * nlocal[1] * NCELL_CUT;
  int proc_zn = get_neighbor(md, 0, 0, -1, off_n);
  int proc_zp = get_neighbor(md, 0, 0, 1, off_p);

  printf("%f %f %f\n", off_n[0], off_n[1], off_n[2]);
  printf("%f %f %f\n", off_p[0], off_p[1], off_p[2]);

  MPI_Irecv(recv_buf_n, ncomm_z, MPI_CHAR, proc_zn, TAG_ZP, mpp->comm, &recv_req_n);
  MPI_Irecv(recv_buf_p, ncomm_z, MPI_CHAR, proc_zp, TAG_ZN, mpp->comm, &recv_req_p);

  esmd_export_box(md, send_buf_n, fields, 0, 0, 0                    , nlocal[0], nlocal[1], NCELL_CUT);
  esmd_export_box(md, send_buf_p, fields, 0, 0, nlocal[2] - NCELL_CUT, nlocal[0], nlocal[1], NCELL_CUT);

  MPI_Isend(send_buf_n, ncomm_z, MPI_CHAR, proc_zn, TAG_ZN, mpp->comm, &send_req_n);
  MPI_Isend(send_buf_p, ncomm_z, MPI_CHAR, proc_zn, TAG_ZP, mpp->comm, &send_req_p);
  MPI_Wait(&recv_req_n, &recv_stat_n);
  MPI_Wait(&recv_req_p, &recv_stat_p);

  esmd_import_box(md, recv_buf_n, fields, 0, 0, -NCELL_CUT, nlocal[0], nlocal[1], NCELL_CUT, off_n);
  esmd_import_box(md, recv_buf_p, fields, 0, 0, nlocal[2] , nlocal[0], nlocal[1], NCELL_CUT, off_p);

  MPI_Wait(&send_req_n, &send_stat_n);
  MPI_Wait(&send_req_p, &send_stat_p);

  
  size_t ncomm_y = entry_size * nlocal[0] * NCELL_CUT * nall[2];
  int proc_yn = get_neighbor(md, 0, -1, 0, off_n);
  int proc_yp = get_neighbor(md, 0, 1, 0, off_p);

  printf("%f %f %f\n", off_n[0], off_n[1], off_n[2]);
  printf("%f %f %f\n", off_p[0], off_p[1], off_p[2]);

  MPI_Irecv(recv_buf_n, ncomm_y, MPI_CHAR, proc_yn, TAG_YP, mpp->comm, &recv_req_n);
  MPI_Irecv(recv_buf_p, ncomm_y, MPI_CHAR, proc_yp, TAG_YN, mpp->comm, &recv_req_p);

  esmd_export_box(md, send_buf_n, fields, 0, 0                    , -NCELL_CUT, nlocal[0], NCELL_CUT, nall[2]);
  esmd_export_box(md, send_buf_p, fields, 0, nlocal[1] - NCELL_CUT, -NCELL_CUT, nlocal[0], NCELL_CUT, nall[2]);

  MPI_Isend(send_buf_n, ncomm_y, MPI_CHAR, proc_yn, TAG_YN, mpp->comm, &send_req_n);
  MPI_Isend(send_buf_p, ncomm_y, MPI_CHAR, proc_yp, TAG_YP, mpp->comm, &send_req_p);

  MPI_Wait(&recv_req_n, &recv_stat_n);
  MPI_Wait(&recv_req_p, &recv_stat_p);

  esmd_import_box(md, recv_buf_n, fields, 0, -NCELL_CUT, -NCELL_CUT, nlocal[0], NCELL_CUT, nall[2], off_n);
  esmd_import_box(md, recv_buf_p, fields, 0, nlocal[1] , -NCELL_CUT, nlocal[0], NCELL_CUT, nall[2], off_p);

  MPI_Wait(&send_req_n, &send_stat_n);
  MPI_Wait(&send_req_p, &send_stat_p);
  
  size_t ncomm_x = entry_size * NCELL_CUT * nall[1] * nall[2];
  int proc_xn = get_neighbor(md, -1, 0, 0, off_n);
  int proc_xp = get_neighbor(md, 1, 0, 0, off_p);

  printf("%f %f %f\n", off_n[0], off_n[1], off_n[2]);
  printf("%f %f %f\n", off_p[0], off_p[1], off_p[2]);

  MPI_Irecv(recv_buf_n, ncomm_x, MPI_CHAR, proc_xn, TAG_XP, mpp->comm, &recv_req_n);
  MPI_Irecv(recv_buf_p, ncomm_x, MPI_CHAR, proc_xp, TAG_XN, mpp->comm, &recv_req_p);

  esmd_export_box(md, send_buf_n, fields, 0                    , -NCELL_CUT, -NCELL_CUT, NCELL_CUT, nall[1], nall[2]);
  esmd_export_box(md, send_buf_p, fields, nlocal[0] - NCELL_CUT, -NCELL_CUT, -NCELL_CUT, NCELL_CUT, nall[1], nall[2]);
  
  MPI_Isend(send_buf_n, ncomm_x, MPI_CHAR, proc_xn, TAG_XN, mpp->comm, &send_req_n);
  MPI_Isend(send_buf_p, ncomm_x, MPI_CHAR, proc_xp, TAG_XP, mpp->comm, &send_req_p);

  MPI_Wait(&recv_req_n, &recv_stat_n);
  MPI_Wait(&recv_req_p, &recv_stat_p);

  esmd_import_box(md, recv_buf_n, fields, -NCELL_CUT, -NCELL_CUT, -NCELL_CUT, NCELL_CUT, nall[1], nall[2], off_n);
  esmd_import_box(md, recv_buf_p, fields, nlocal[0] , -NCELL_CUT, -NCELL_CUT, NCELL_CUT, nall[1], nall[2], off_p);

  MPI_Wait(&send_req_n, &send_stat_n);
  MPI_Wait(&send_req_p, &send_stat_p);


  esmd_free(send_buf_n);
  esmd_free(send_buf_p);
  esmd_free(recv_buf_n);
  esmd_free(recv_buf_p);
}
