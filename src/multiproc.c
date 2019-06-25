#include <mpi.h>
#include <data.h>
#include <util.h>
#include <stdio.h>
#include <assert.h>

#include <multiproc.h>
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

static inline void init_halo_ordered(halo_t *halo, esmd_t *md, int dir, int del, void *send_buf, void *recv_buf) {
  for (int i = 0; i < dir; i ++) {
    halo->off[i][0] = 0;
    halo->off[i][1] = 0;
    halo->len[i] = md->box.nlocal[i];
  }

  for (int i = dir + 1; i < 3; i ++) {
    halo->off[i][0] = -NCELL_CUT;
    halo->off[i][1] = -NCELL_CUT;
    halo->len[i] = md->box.nall[i];
  }
  halo->off[dir][0] = del < 0 ? 0 : md->box.nlocal[dir] - NCELL_CUT;
  halo->off[dir][1] = del < 0 ? -NCELL_CUT : md->box.nlocal[dir];
  halo->len[dir] = NCELL_CUT;
  halo->ncells = halo->len[0] * halo->len[1] * halo->len[2];
  int del3[3] = {0, 0, 0};
  del3[dir] = del;
  halo->neighbor = get_neighbor(md, del3[0], del3[1], del3[2], halo->translation);
  halo->send_buf = send_buf;
  halo->recv_buf = recv_buf;
  halo->send_tag = 100 + dir * 10 + del;
  halo->recv_tag = 100 + dir * 10 - del;
}

void init_comm_ordered(esmd_t *md){
  int *nall = md->box.nall;
  size_t max_entry_size = sizeof(cell_t) + sizeof(celldata_t);
  int max_comm_cells = max(nall[0] * nall[1], max(nall[1] * nall[2], nall[0] * nall[2])) * NCELL_CUT;

  void *send_buf_next = esmd_malloc(max_comm_cells * max_entry_size, "exchange_buffer");
  void *send_buf_prev = esmd_malloc(max_comm_cells * max_entry_size, "exchange_buffer");
  void *recv_buf_next = esmd_malloc(max_comm_cells * max_entry_size, "exchange_buffer");
  void *recv_buf_prev = esmd_malloc(max_comm_cells * max_entry_size, "exchange_buffer");

  halo_t *halo = md->mpp.halo;

  init_halo_ordered(halo + 0, md, 0, -1, send_buf_prev, recv_buf_prev);      
  init_halo_ordered(halo + 1, md, 0,  1, send_buf_next, recv_buf_next);
  init_halo_ordered(halo + 2, md, 1, -1, send_buf_prev, recv_buf_prev);
  init_halo_ordered(halo + 3, md, 1,  1, send_buf_next, recv_buf_next);
  init_halo_ordered(halo + 4, md, 2, -1, send_buf_prev, recv_buf_prev);
  init_halo_ordered(halo + 5, md, 2,  1, send_buf_next, recv_buf_next);
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
  
  box->nall[0] = cntx + NCELL_CUT * 2;
  box->nall[1] = cnty + NCELL_CUT * 2;
  box->nall[2] = cntz + NCELL_CUT * 2;
  
  box->nlocal[0] = cntx;
  box->nlocal[1] = cnty;
  box->nlocal[2] = cntz;

  box->offset[0] = stx;
  box->offset[1] = sty;
  box->offset[2] = stz;

  init_comm_ordered(md);
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

void esmd_exchange_cell_local_to_halo(esmd_t *md, int fields, int flags) {
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

  MPI_Irecv(recv_buf_n, ncomm_z, MPI_CHAR, proc_zn, TAG_ZP, mpp->comm, &recv_req_n);
  MPI_Irecv(recv_buf_p, ncomm_z, MPI_CHAR, proc_zp, TAG_ZN, mpp->comm, &recv_req_p);

  esmd_export_box(md, send_buf_n, fields, flags, 0, 0, 0                    , nlocal[0], nlocal[1], NCELL_CUT);
  esmd_export_box(md, send_buf_p, fields, flags, 0, 0, nlocal[2] - NCELL_CUT, nlocal[0], nlocal[1], NCELL_CUT);

  MPI_Isend(send_buf_n, ncomm_z, MPI_CHAR, proc_zn, TAG_ZN, mpp->comm, &send_req_n);
  MPI_Isend(send_buf_p, ncomm_z, MPI_CHAR, proc_zn, TAG_ZP, mpp->comm, &send_req_p);
  MPI_Wait(&recv_req_n, &recv_stat_n);
  MPI_Wait(&recv_req_p, &recv_stat_p);

  esmd_import_box(md, recv_buf_n, fields, flags, 0, 0, -NCELL_CUT, nlocal[0], nlocal[1], NCELL_CUT, off_n);
  esmd_import_box(md, recv_buf_p, fields, flags, 0, 0, nlocal[2] , nlocal[0], nlocal[1], NCELL_CUT, off_p);

  MPI_Wait(&send_req_n, &send_stat_n);
  MPI_Wait(&send_req_p, &send_stat_p);

  
  size_t ncomm_y = entry_size * nlocal[0] * NCELL_CUT * nall[2];
  int proc_yn = get_neighbor(md, 0, -1, 0, off_n);
  int proc_yp = get_neighbor(md, 0, 1, 0, off_p);

  MPI_Irecv(recv_buf_n, ncomm_y, MPI_CHAR, proc_yn, TAG_YP, mpp->comm, &recv_req_n);
  MPI_Irecv(recv_buf_p, ncomm_y, MPI_CHAR, proc_yp, TAG_YN, mpp->comm, &recv_req_p);

  esmd_export_box(md, send_buf_n, fields, flags, 0, 0                    , -NCELL_CUT, nlocal[0], NCELL_CUT, nall[2]);
  esmd_export_box(md, send_buf_p, fields, flags, 0, nlocal[1] - NCELL_CUT, -NCELL_CUT, nlocal[0], NCELL_CUT, nall[2]);

  MPI_Isend(send_buf_n, ncomm_y, MPI_CHAR, proc_yn, TAG_YN, mpp->comm, &send_req_n);
  MPI_Isend(send_buf_p, ncomm_y, MPI_CHAR, proc_yp, TAG_YP, mpp->comm, &send_req_p);

  MPI_Wait(&recv_req_n, &recv_stat_n);
  MPI_Wait(&recv_req_p, &recv_stat_p);

  esmd_import_box(md, recv_buf_n, fields, flags, 0, -NCELL_CUT, -NCELL_CUT, nlocal[0], NCELL_CUT, nall[2], off_n);
  esmd_import_box(md, recv_buf_p, fields, flags, 0, nlocal[1] , -NCELL_CUT, nlocal[0], NCELL_CUT, nall[2], off_p);

  MPI_Wait(&send_req_n, &send_stat_n);
  MPI_Wait(&send_req_p, &send_stat_p);
  
  size_t ncomm_x = entry_size * NCELL_CUT * nall[1] * nall[2];
  int proc_xn = get_neighbor(md, -1, 0, 0, off_n);
  int proc_xp = get_neighbor(md, 1, 0, 0, off_p);

  printf("%f %f %f\n", off_n[0], off_n[1], off_n[2]);
  printf("%f %f %f\n", off_p[0], off_p[1], off_p[2]);

  MPI_Irecv(recv_buf_n, ncomm_x, MPI_CHAR, proc_xn, TAG_XP, mpp->comm, &recv_req_n);
  MPI_Irecv(recv_buf_p, ncomm_x, MPI_CHAR, proc_xp, TAG_XN, mpp->comm, &recv_req_p);

  esmd_export_box(md, send_buf_n, fields, flags, 0                    , -NCELL_CUT, -NCELL_CUT, NCELL_CUT, nall[1], nall[2]);
  esmd_export_box(md, send_buf_p, fields, flags, nlocal[0] - NCELL_CUT, -NCELL_CUT, -NCELL_CUT, NCELL_CUT, nall[1], nall[2]);
  
  MPI_Isend(send_buf_n, ncomm_x, MPI_CHAR, proc_xn, TAG_XN, mpp->comm, &send_req_n);
  MPI_Isend(send_buf_p, ncomm_x, MPI_CHAR, proc_xp, TAG_XP, mpp->comm, &send_req_p);

  MPI_Wait(&recv_req_n, &recv_stat_n);
  MPI_Wait(&recv_req_p, &recv_stat_p);

  esmd_import_box(md, recv_buf_n, fields, flags, -NCELL_CUT, -NCELL_CUT, -NCELL_CUT, NCELL_CUT, nall[1], nall[2], off_n);
  esmd_import_box(md, recv_buf_p, fields, flags, nlocal[0] , -NCELL_CUT, -NCELL_CUT, NCELL_CUT, nall[1], nall[2], off_p);

  MPI_Wait(&send_req_n, &send_stat_n);
  MPI_Wait(&send_req_p, &send_stat_p);

  esmd_free(send_buf_n);
  esmd_free(send_buf_p);
  esmd_free(recv_buf_n);
  esmd_free(recv_buf_p);
}

void esmd_exchange_cell_halo_to_local(esmd_t *md, int fields, int flags) {
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

  size_t ncomm_x = entry_size * NCELL_CUT * nall[1] * nall[2];
  int proc_xn = get_neighbor(md, -1, 0, 0, off_n);
  int proc_xp = get_neighbor(md, 1, 0, 0, off_p);

  MPI_Irecv(recv_buf_n, ncomm_x, MPI_CHAR, proc_xn, TAG_XP, mpp->comm, &recv_req_n);
  MPI_Irecv(recv_buf_p, ncomm_x, MPI_CHAR, proc_xp, TAG_XN, mpp->comm, &recv_req_p);

  esmd_export_box(md, send_buf_n, fields, flags, -NCELL_CUT, -NCELL_CUT, -NCELL_CUT, NCELL_CUT, nall[1], nall[2]);
  esmd_export_box(md, send_buf_p, fields, flags, nlocal[0] , -NCELL_CUT, -NCELL_CUT, NCELL_CUT, nall[1], nall[2]);
  
  MPI_Isend(send_buf_n, ncomm_x, MPI_CHAR, proc_xn, TAG_XN, mpp->comm, &send_req_n);
  MPI_Isend(send_buf_p, ncomm_x, MPI_CHAR, proc_xp, TAG_XP, mpp->comm, &send_req_p);

  MPI_Wait(&recv_req_n, &recv_stat_n);
  MPI_Wait(&recv_req_p, &recv_stat_p);

  esmd_import_box(md, recv_buf_n, fields, flags, 0                    , -NCELL_CUT, -NCELL_CUT, NCELL_CUT, nall[1], nall[2], off_n);
  esmd_import_box(md, recv_buf_p, fields, flags, nlocal[0] - NCELL_CUT, -NCELL_CUT, -NCELL_CUT, NCELL_CUT, nall[1], nall[2], off_p);

  MPI_Wait(&send_req_n, &send_stat_n);
  MPI_Wait(&send_req_p, &send_stat_p);

  size_t ncomm_y = entry_size * nlocal[0] * NCELL_CUT * nall[2];
  int proc_yn = get_neighbor(md, 0, -1, 0, off_n);
  int proc_yp = get_neighbor(md, 0, 1, 0, off_p);

  MPI_Irecv(recv_buf_n, ncomm_y, MPI_CHAR, proc_yn, TAG_YP, mpp->comm, &recv_req_n);
  MPI_Irecv(recv_buf_p, ncomm_y, MPI_CHAR, proc_yp, TAG_YN, mpp->comm, &recv_req_p);

  esmd_export_box(md, send_buf_n, fields, flags, 0, -NCELL_CUT, -NCELL_CUT, nlocal[0], NCELL_CUT, nall[2]);
  esmd_export_box(md, send_buf_p, fields, flags, 0, nlocal[1] , -NCELL_CUT, nlocal[0], NCELL_CUT, nall[2]);

  MPI_Isend(send_buf_n, ncomm_y, MPI_CHAR, proc_yn, TAG_YN, mpp->comm, &send_req_n);
  MPI_Isend(send_buf_p, ncomm_y, MPI_CHAR, proc_yp, TAG_YP, mpp->comm, &send_req_p);

  MPI_Wait(&recv_req_n, &recv_stat_n);
  MPI_Wait(&recv_req_p, &recv_stat_p);

  esmd_import_box(md, recv_buf_n, fields, flags, 0, 0                    , -NCELL_CUT, nlocal[0], NCELL_CUT, nall[2], off_n);
  esmd_import_box(md, recv_buf_p, fields, flags, 0, nlocal[1] - NCELL_CUT, -NCELL_CUT, nlocal[0], NCELL_CUT, nall[2], off_p);

  MPI_Wait(&send_req_n, &send_stat_n);
  MPI_Wait(&send_req_p, &send_stat_p);

  size_t ncomm_z = entry_size * nlocal[0] * nlocal[1] * NCELL_CUT;
  int proc_zn = get_neighbor(md, 0, 0, -1, off_n);
  int proc_zp = get_neighbor(md, 0, 0, 1, off_p);

  MPI_Irecv(recv_buf_n, ncomm_z, MPI_CHAR, proc_zn, TAG_ZP, mpp->comm, &recv_req_n);
  MPI_Irecv(recv_buf_p, ncomm_z, MPI_CHAR, proc_zp, TAG_ZN, mpp->comm, &recv_req_p);

  esmd_export_box(md, send_buf_n, fields, flags, 0, 0, -NCELL_CUT, nlocal[0], nlocal[1], NCELL_CUT);
  esmd_export_box(md, send_buf_p, fields, flags, 0, 0, nlocal[2] , nlocal[0], nlocal[1], NCELL_CUT);

  MPI_Isend(send_buf_n, ncomm_z, MPI_CHAR, proc_zn, TAG_ZN, mpp->comm, &send_req_n);
  MPI_Isend(send_buf_p, ncomm_z, MPI_CHAR, proc_zn, TAG_ZP, mpp->comm, &send_req_p);
  MPI_Wait(&recv_req_n, &recv_stat_n);
  MPI_Wait(&recv_req_p, &recv_stat_p);

  esmd_import_box(md, recv_buf_n, fields, flags, 0, 0, 0                    , nlocal[0], nlocal[1], NCELL_CUT, off_n);
  esmd_import_box(md, recv_buf_p, fields, flags, 0, 0, nlocal[2] - NCELL_CUT, nlocal[0], nlocal[1], NCELL_CUT, off_p);

  MPI_Wait(&send_req_n, &send_stat_n);
  MPI_Wait(&send_req_p, &send_stat_p);  


  esmd_free(send_buf_n);
  esmd_free(send_buf_p);
  esmd_free(recv_buf_n);
  esmd_free(recv_buf_p);
}

void esmd_comm_start(esmd_t *md, multiproc_t *mpp, halo_t *halo, int dir, int fields, int flags){
  size_t entry_size = esmd_fields_size(fields);
  int neigh = halo->neighbor;
  size_t comm_size = halo->ncells * entry_size;
  
  int off_x = halo->off[0][dir], off_y = halo->off[1][dir], off_z = halo->off[2][dir];
  int len_x = halo->len[0], len_y = halo->len[1], len_z = halo->len[2];
  printf("send to %d, box={%d, %d, %d, %d, %d, %d}\n", neigh, off_x, off_y, off_z, len_x, len_y, len_z);
  MPI_Irecv(halo->recv_buf, comm_size, MPI_CHAR, neigh, halo->recv_tag, mpp->comm, &halo->recv_req);
  esmd_export_box(md, halo->send_buf, fields, flags, off_x, off_y, off_z, len_x, len_y, len_z);
  MPI_Isend(halo->send_buf, comm_size, MPI_CHAR, neigh, halo->send_tag, mpp->comm, &halo->send_req);
}

void esmd_comm_finish(esmd_t *md, multiproc_t *mpp, halo_t *halo, int dir, int fields, int flags){
  int off_x = halo->off[0][1 - dir], off_y = halo->off[1][1 - dir], off_z = halo->off[2][1 - dir];
  int len_x = halo->len[0], len_y = halo->len[1], len_z = halo->len[2];
  areal *trans = halo->translation;
  printf("recv from %d, box={%d, %d, %d, %d, %d, %d}\n", halo->neighbor, off_x, off_y, off_z, len_x, len_y, len_z);
  MPI_Wait(&halo->recv_req, &halo->recv_stat);
  esmd_import_box(md, halo->recv_buf, fields, flags, off_x, off_y, off_z, len_x, len_y, len_z, trans);
  MPI_Wait(&halo->send_req, &halo->send_stat);
}

void esmd_exchange_cell(esmd_t *md, int direction, int fields, int flags) {
  box_t *box = &(md->box);
  multiproc_t *mpp = &(md->mpp);
  if (direction == LOCAL_TO_HALO) {
    esmd_comm_start(md, mpp, mpp->halo + 4, direction, fields, flags);
    esmd_comm_start(md, mpp, mpp->halo + 5, direction, fields, flags);
    esmd_comm_finish(md, mpp, mpp->halo + 4, direction, fields, flags);
    esmd_comm_finish(md, mpp, mpp->halo + 5, direction, fields, flags);

    esmd_comm_start(md, mpp, mpp->halo + 2, direction, fields, flags);
    esmd_comm_start(md, mpp, mpp->halo + 3, direction, fields, flags);
    esmd_comm_finish(md, mpp, mpp->halo + 2, direction, fields, flags);
    esmd_comm_finish(md, mpp, mpp->halo + 3, direction, fields, flags);
    
    esmd_comm_start(md, mpp, mpp->halo + 0, direction, fields, flags);
    esmd_comm_start(md, mpp, mpp->halo + 1, direction, fields, flags);
    esmd_comm_finish(md, mpp, mpp->halo + 0, direction, fields, flags);
    esmd_comm_finish(md, mpp, mpp->halo + 1, direction, fields, flags);
  } else if (direction == HALO_TO_LOCAL){
    esmd_comm_start(md, mpp, mpp->halo + 0, direction, fields, flags);
    esmd_comm_start(md, mpp, mpp->halo + 1, direction, fields, flags);
    esmd_comm_finish(md, mpp, mpp->halo + 0, direction, fields, flags);
    esmd_comm_finish(md, mpp, mpp->halo + 1, direction, fields, flags);

    esmd_comm_start(md, mpp, mpp->halo + 2, direction, fields, flags);
    esmd_comm_start(md, mpp, mpp->halo + 3, direction, fields, flags);
    esmd_comm_finish(md, mpp, mpp->halo + 2, direction, fields, flags);
    esmd_comm_finish(md, mpp, mpp->halo + 3, direction, fields, flags);

    esmd_comm_start(md, mpp, mpp->halo + 4, direction, fields, flags);
    esmd_comm_start(md, mpp, mpp->halo + 5, direction, fields, flags);
    esmd_comm_finish(md, mpp, mpp->halo + 4, direction, fields, flags);
    esmd_comm_finish(md, mpp, mpp->halo + 5, direction, fields, flags);
  }
}


