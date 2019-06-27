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

void init_halo_unordered(halo_t *halo, esmd_t *md, int dx, int dy, int dz) {
  int d[3] = {dx, dy, dz};
  for (int i = 0; i < 3; i ++) {
    if (d[i] == 0){
      halo->off[i][0] = 0;
      halo->off[i][1] = 0;
      halo->len[i] = md->box.nlocal[i];
    } else {
      halo->off[i][0] = (d[i] == -1) ? 0 : md->box.nlocal[i] - NCELL_CUT;
      halo->off[i][1] = (d[i] == -1) ? -NCELL_CUT : md->box.nlocal[i];
      halo->len[i] = NCELL_CUT;
    }
  }
  halo->ncells = halo->len[0] * halo->len[1] * halo->len[2];
  halo->neighbor = get_neighbor(md, dx, dy, dz, halo->translation);
  //printf("%d %d %d %d %d %d %d %d\n", dx, dy, dz, halo->neighbor, halo->ncells, halo->len[0], halo->len[1], halo->len[2]);
  size_t max_entry_size = sizeof(cell_t) + sizeof(celldata_t);
  halo->send_buf = esmd_malloc(halo->ncells * max_entry_size, "exchange_buffer");
  halo->recv_buf = esmd_malloc(halo->ncells * max_entry_size, "exchange_buffer");
  halo->send_tag = 200 + dx * 9 + dy * 3 + dz;
  halo->recv_tag = 200 - dx * 9 - dy * 3 - dz ;
    
}

void init_comm_unordered(esmd_t *md){
  int nhalo = 0;
  for (int dx = -1; dx <= 1; dx ++) {
    for (int dy = -1; dy <= 1; dy ++) {
      for (int dz = -1; dz <= 1; dz ++) {
        if (dx == 0 && dy == 0 && dz == 0) continue;
        init_halo_unordered(md->mpp.halo + nhalo, md, dx, dy, dz);
        nhalo ++;
      }
    }
  }
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

  init_comm_unordered(md);
}

inline int proc_3d_to_flat(multiproc_t *mpp, int pidx, int pidy, int pidz){
  return pidx * mpp->npy * mpp->npz + pidy * mpp->npz + pidz;
}


void esmd_comm_start(esmd_t *md, MPI_Comm comm, halo_t *halo, int dir, int fields, int flags){
  size_t entry_size = esmd_fields_size(fields);
  int neigh = halo->neighbor;
  size_t comm_size = halo->ncells * entry_size;
  
  int off_x = halo->off[0][dir], off_y = halo->off[1][dir], off_z = halo->off[2][dir];
  int len_x = halo->len[0], len_y = halo->len[1], len_z = halo->len[2];

  MPI_Irecv(halo->recv_buf, comm_size, MPI_CHAR, neigh, halo->recv_tag, comm, &halo->recv_req);
  esmd_export_box(md, halo->send_buf, fields, flags, off_x, off_y, off_z, len_x, len_y, len_z);
  MPI_Isend(halo->send_buf, comm_size, MPI_CHAR, neigh, halo->send_tag, comm, &halo->send_req);
}

void esmd_comm_finish(esmd_t *md, halo_t *halo, int dir, int fields, int flags){
  int off_x = halo->off[0][1 - dir], off_y = halo->off[1][1 - dir], off_z = halo->off[2][1 - dir];
  int len_x = halo->len[0], len_y = halo->len[1], len_z = halo->len[2];
  areal *trans = halo->translation;

  MPI_Wait(&halo->recv_req, &halo->recv_stat);
  esmd_import_box(md, halo->recv_buf, fields, flags, off_x, off_y, off_z, len_x, len_y, len_z, trans);
  MPI_Wait(&halo->send_req, &halo->send_stat);
}

void esmd_exchange_cell(esmd_t *md, int direction, int fields, int flags) {
  for (int i = 0; i < 26; i ++){
    //printf("%d\n", md->mpp.halo[i].neighbor);
    esmd_comm_start(md, md->mpp.comm, md->mpp.halo + i, direction, fields, flags);
  }
  for (int i = 0; i < 26; i ++){
    esmd_comm_finish(md, md->mpp.halo + i, direction, fields, flags);
  }
}
void esmd_exchange_cell_ordered(esmd_t *md, int direction, int fields, int flags) {
  multiproc_t *mpp = &(md->mpp);
  if (direction == LOCAL_TO_HALO) {
    for (int i = 4; i >= 0; i -= 2){
      //printf("%d %d\n", mpp->halo[i + 0].neighbor, mpp->halo[i + 1].neighbor);
      esmd_comm_start(md, mpp->comm, mpp->halo + i + 0, direction, fields, flags);
      esmd_comm_start(md, mpp->comm, mpp->halo + i + 1, direction, fields, flags);
      esmd_comm_finish(md, mpp->halo + i + 0, direction, fields, flags);
      esmd_comm_finish(md, mpp->halo + i + 1, direction, fields, flags);
      //printf("%d %d\n", mpp->halo[i + 0].neighbor, mpp->halo[i + 1].neighbor);
    }
  } else if (direction == HALO_TO_LOCAL){
    for (int i = 0; i < 6; i += 2){
      //printf("%d %d\n", mpp->halo[i + 0].neighbor, mpp->halo[i + 1].neighbor);
      esmd_comm_start(md, mpp->comm, mpp->halo + i + 0, direction, fields, flags);
      esmd_comm_start(md, mpp->comm, mpp->halo + i + 1, direction, fields, flags);
      esmd_comm_finish(md, mpp->halo + i + 0, direction, fields, flags);
      esmd_comm_finish(md, mpp->halo + i + 1, direction, fields, flags);
      //printf("%d %d\n", mpp->halo[i + 0].neighbor, mpp->halo[i + 1].neighbor);
    }
  }
  
}
