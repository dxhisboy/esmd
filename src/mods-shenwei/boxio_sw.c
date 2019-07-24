#include <data.h>
#include <util.h>
#define DEBUG_THIS_FILE
#include <log.h>
#define NPES_IO 64
typedef struct boxio_param {
  box_t *box;
  void *buffer;
  int offset_atoms[NPES_IO + 1];
  int offset_cells[NPES_IO + 1];
  int flags, fields;
  int xlo, ylo, zlo;
  int xlen, ylen, zlen;
} boxio_param_t;
inline size_t estimate_atoms_size(int fields){
  size_t accu_size = 0;
  if (fields & CELL_X) accu_size += sizeof(areal) * 3;
  if (fields & CELL_V) accu_size += sizeof(areal) * 3;
  if (fields & CELL_F) accu_size += sizeof(areal) * 3;
  if (fields & CELL_Q) accu_size += sizeof(ireal);
  if (fields & CELL_T) accu_size += sizeof(int);
  if (fields & CELL_E) accu_size += sizeof(int);
  return accu_size;
}
#ifdef CPE
#include <slave.h>
#include <dma_macros.h>
#define BUF_SIZE 2048
#define export_void(start, size) {                      \
    if (l_offset + (size) > BUF_SIZE) {                 \
      pe_put(g_buf + g_offset, l_buf, l_offset);        \
      dma_syn();                                        \
      g_offset += l_offset;                             \
      l_offset = 0;                                     \
    }                                                   \
    pe_get((start), l_buf + l_offset, (size));          \
    dma_syn();                                          \
    l_offset += (size);                                 \
  }

void export_box_cpe(boxio_param_t *gl_pm){
  dma_init();
  boxio_param_t pm;
  pe_get(gl_pm, &pm, sizeof(boxio_param_t));
  dma_syn();
  box_t box;
  pe_get(pm.box, &box, sizeof(box_t));
  dma_syn();
  char l_buf[BUF_SIZE];
  size_t l_offset = 0;
  size_t g_offset = 0;

  int xlen = pm.xlen, ylen = pm.ylen, zlen = pm.zlen;
  int xlo = pm.xlo, ylo = pm.ylo, zlo = pm.zlo;

  void *g_buf = pm.buffer;
  int flags = pm.flags;
  int fields = pm.fields;
  int *offset_cells = pm.offset_cells;

  if (fields & CELL_META) g_offset += offset_cells[_MYID] * sizeof(cell_t);
  size_t atomdata_size = estimate_atoms_size(pm.fields);
  g_offset += atomdata_size * pm.offset_atoms[_MYID];

  for (int icell = offset_cells[_MYID]; icell < offset_cells[_MYID + 1]; icell ++){
    int ii = icell % xlen + xlo;
    int jj = icell / xlen % ylen + ylo;
    int kk = icell / xlen / ylen + zlo;
    int celloff = get_cell_off((&box), ii, jj, kk);
    cell_t cell;

    if (fields & CELL_META){
      export_void(box.cells + celloff, sizeof(cell_t));
      memcpy(&cell, l_buf + l_offset - sizeof(cell_t), sizeof(cell_t));
    } else {
      pe_get(box.cells + celloff, &cell, sizeof(cell_t));
      dma_syn();
    }
    celldata_t *data = box.celldata + celloff;
    if ((flags & TRANS_ATOMS) && cell.natoms > 0){
      int off = 0;
      int cnt = cell.natoms;
      if (fields & CELL_X) export_void(data->x + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_V) export_void(data->v + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_F) export_void(data->f + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_Q) export_void(data->q + off, sizeof(ireal) * cnt);
      if (fields & CELL_T) export_void(data->type + off, sizeof(int) * cnt);
      if (fields & CELL_E) export_void(data->export + off, sizeof(int) * cnt);
    }
    if ((flags & TRANS_EXPORTS) && cell.export_ptr != CELL_SIZE){
      int off = cell.export_ptr;
      int cnt = CELL_SIZE - off;
      if (fields & CELL_X) export_void(data->x + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_V) export_void(data->v + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_F) export_void(data->f + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_Q) export_void(data->q + off, sizeof(ireal) * cnt);
      if (fields & CELL_T) export_void(data->type + off, sizeof(int) * cnt);
      if (fields & CELL_E) export_void(data->export + off, sizeof(int) * cnt);
    }
    if (!(flags & (TRANS_ATOMS | TRANS_EXPORTS))){
      int off = 0;
      int cnt = CELL_SIZE;
      if (fields & CELL_X) export_void(data->x + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_V) export_void(data->v + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_F) export_void(data->f + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_Q) export_void(data->q + off, sizeof(ireal) * cnt);
      if (fields & CELL_T) export_void(data->type + off, sizeof(int) * cnt);
      if (fields & CELL_E) export_void(data->export + off, sizeof(int) * cnt);
    }
  }
  if (l_offset) {
    pe_put(g_buf + g_offset, l_buf, l_offset);
    dma_syn();
  }
}

#define import_void(start, size) {                                      \
    if (l_offset + (size) > BUF_SIZE) {                                 \
      memcpy(l_buf, l_buf + l_offset, BUF_SIZE - l_offset);             \
      pe_get(g_buf + g_offset, l_buf + BUF_SIZE - l_offset, l_offset);  \
      dma_syn();                                                        \
      g_offset += l_offset;                                             \
      l_offset = 0;                                                     \
    }                                                                   \
    pe_put((start), l_buf + l_offset, (size));                          \
    dma_syn();                                                          \
    l_offset += (size);                                                 \
  }
#undef BUF_SIZE
#define BUF_SIZE 2048
void import_box_cpe(boxio_param_t *gl_pm){
  dma_init();
  boxio_param_t pm;
  pe_get(gl_pm, &pm, sizeof(boxio_param_t));
  dma_syn();
  box_t box;
  pe_get(pm.box, &box, sizeof(box_t));
  dma_syn();
  char l_buf[BUF_SIZE];
  size_t l_offset = BUF_SIZE;
  size_t g_offset = 0;

  int xlen = pm.xlen, ylen = pm.ylen, zlen = pm.zlen;
  int xlo = pm.xlo, ylo = pm.ylo, zlo = pm.zlo;

  void *g_buf = pm.buffer;
  int flags = pm.flags;
  int fields = pm.fields;
  int *offset_cells = pm.offset_cells;

  if (fields & CELL_META) g_offset += offset_cells[_MYID] * sizeof(cell_t);
  size_t atomdata_size = estimate_atoms_size(pm.fields);
  g_offset += atomdata_size * pm.offset_atoms[_MYID];
  //if (_MYID != 8) return;
  for (int icell = offset_cells[_MYID]; icell < offset_cells[_MYID + 1]; icell ++){
    int ii = icell % xlen + xlo;
    int jj = icell / xlen % ylen + ylo;
    int kk = icell / xlen / ylen + zlo;
    int celloff = get_cell_off((&box), ii, jj, kk);
    cell_t cell;

    if (fields & CELL_META){
      import_void(box.cells + celloff, sizeof(cell_t));
      memcpy(&cell, l_buf + l_offset - sizeof(cell_t), sizeof(cell_t));
    } else {
      pe_get(box.cells + celloff, &cell, sizeof(cell_t));
      dma_syn();
    }
    celldata_t *data = box.celldata + celloff;
    if ((flags & TRANS_ATOMS) && cell.natoms > 0){
      int off = 0;
      int cnt = cell.natoms;
      //printf("%p %d\n", data->x, cnt);
      if (fields & CELL_X) import_void(data->x + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_V) import_void(data->v + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_F) import_void(data->f + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_Q) import_void(data->q + off, sizeof(ireal) * cnt);
      if (fields & CELL_T) import_void(data->type + off, sizeof(int) * cnt);
      if (fields & CELL_E) import_void(data->export + off, sizeof(int) * cnt);
    }
    if ((flags & TRANS_EXPORTS) && cell.export_ptr != CELL_SIZE){
      int off = cell.export_ptr;
      int cnt = CELL_SIZE - off;
      if (fields & CELL_X) import_void(data->x + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_V) import_void(data->v + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_F) import_void(data->f + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_Q) import_void(data->q + off, sizeof(ireal) * cnt);
      if (fields & CELL_T) import_void(data->type + off, sizeof(int) * cnt);
      if (fields & CELL_E) import_void(data->export + off, sizeof(int) * cnt);
    }
    if (!(flags & (TRANS_ATOMS | TRANS_EXPORTS))){
      int off = 0;
      int cnt = CELL_SIZE;
      if (fields & CELL_X) import_void(data->x + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_V) import_void(data->v + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_F) import_void(data->f + off, sizeof(areal) * 3 * cnt);
      if (fields & CELL_Q) import_void(data->q + off, sizeof(ireal) * cnt);
      if (fields & CELL_T) import_void(data->type + off, sizeof(int) * cnt);
      if (fields & CELL_E) import_void(data->export + off, sizeof(int) * cnt);
    }
  }
  /* if (l_offset) { */
  /*   pe_put(g_buf + g_offset, l_buf, l_offset); */
  /*   dma_syn(); */
  /* } */
}
#endif

#ifdef MPE
#include <athread.h>
#include <timer.h>
extern void slave_export_box_cpe(boxio_param_t *);
extern void slave_import_box_cpe(boxio_param_t *);
inline void
estimate_io_size(esmd_t *md, int *offset_atoms, int *offset_cells, int flags,
                 int xlo, int ylo, int zlo, int xlen, int ylen, int zlen){
  box_t *box = md->box;
  for (int pe = 0; pe <= NPES_IO; pe ++){
    offset_atoms[pe] = 0;
    offset_cells[pe] = 0;
  }
  int ncells = xlen * ylen * zlen;
  int ncell_div = ncells / NPES_IO;
  int ncell_rem = ncells % NPES_IO;
  for (int pe = 0; pe < NPES_IO; pe ++){
    offset_cells[pe + 1] = ncell_div;
    if (pe < ncell_rem)
      offset_cells[pe + 1] ++;
    offset_cells[pe + 1] += offset_cells[pe];
  }
  for (int pe = 0; pe < NPES_IO; pe ++){
    for (int icell = offset_cells[pe]; icell < offset_cells[pe + 1]; icell ++){
      int ii = icell % xlen + xlo;
      int jj = icell / xlen % ylen + ylo;
      int kk = icell / (xlen * ylen) + zlo;
      int celloff = get_cell_off(box, ii, jj, kk);
      if (abs(box->cells[celloff].natoms) > 100)
        debug("%d %d %d %d %d\n", ii, jj, kk, icell, ncells);
      if (flags & (TRANS_ATOMS | TRANS_EXPORTS)){
        if (flags & TRANS_ATOMS)
          offset_atoms[pe + 1] += box->cells[celloff].natoms;
        if (flags & TRANS_EXPORTS)
          offset_atoms[pe + 1] += CELL_SIZE - box->cells[celloff].export_ptr;
      } else {
        offset_atoms[pe + 1] += CELL_SIZE;
      }
    }
    offset_atoms[pe + 1] += offset_atoms[pe];
  }
  /* for (int i = 0; i <= 64; i ++){ */
  /*   debug("%d %d %d\n", ncells, offset_atoms[i], offset_cells[i]); */
  /* } */
}

size_t esmd_export_box_sw(esmd_t *md, void *buffer, int fields, int flags, int xlo, int ylo, int zlo, int xlen, int ylen, int zlen){
  timer_start("export_box_sw");
  boxio_param_t pm;
  estimate_io_size(md, pm.offset_atoms, pm.offset_cells, flags, xlo, ylo, zlo, xlen, ylen, zlen);
  pm.xlen = xlen;
  pm.ylen = ylen;
  pm.zlen = zlen;
  pm.xlo = xlo;
  pm.ylo = ylo;
  pm.zlo = zlo;
  pm.flags = flags;
  pm.fields = fields;
  pm.buffer = buffer + sizeof(int) * (NPES_IO + 1) * 2;
  pm.box = md->box;
  athread_spawn(export_box_cpe, &pm);
  memcpy(buffer, pm.offset_atoms, sizeof(int) * (NPES_IO + 1));
  memcpy(buffer + sizeof(int) * (NPES_IO + 1), pm.offset_cells, sizeof(int) * (NPES_IO + 1));
  size_t ret = pm.offset_atoms[NPES_IO] * estimate_atoms_size(fields);
  if (fields & CELL_META) ret += pm.offset_cells[NPES_IO] * sizeof(cell_t);
  ret += sizeof(int) * (NPES_IO + 1) * 2;
  athread_join();
  timer_stop("export_box_sw");
  //printf("%ld %ld %ld %ld\n", ret, pm.offset_atoms[NPES_IO], pm.offset_cells[NPES_IO], estimate_atoms_size(fields));
  return ret;
}

size_t esmd_import_box_sw(esmd_t *md, void *buffer, int fields, int flags, int xlo, int ylo, int zlo, int xlen, int ylen, int zlen, areal *off){
  timer_start("import_box_sw");
  void *buffer_head = buffer + sizeof(int) * (NPES_IO + 1) * 2;
  //esmd_import_box(md, buffer_head, fields, flags, xlo, ylo, zlo, xlen, ylen, zlen, off);
  boxio_param_t pm;
  //estimate_io_size(md, pm.offset_atoms, pm.offset_cells, flags, xlo, ylo, zlo, xlen, ylen, zlen);
  memcpy(pm.offset_atoms, buffer, sizeof(int) * (NPES_IO + 1));
  memcpy(pm.offset_cells, buffer + sizeof(int) * (NPES_IO + 1), sizeof(int) * (NPES_IO + 1));
  pm.xlen = xlen;
  pm.ylen = ylen;
  pm.zlen = zlen;
  pm.xlo = xlo;
  pm.ylo = ylo;
  pm.zlo = zlo;
  pm.flags = flags;
  pm.fields = fields;
  pm.buffer = buffer_head;
  pm.box = md->box;
  /* for (int i = 0; i <= 64; i ++){ */
  /*   debug("%d %d %d\n", xlen * ylen * zlen, pm.offset_atoms[i], pm.offset_cells[i]); */
  /* } */

  athread_spawn(import_box_cpe, &pm);
  size_t ret = pm.offset_atoms[NPES_IO] * estimate_atoms_size(fields);
  if (fields & CELL_META) ret += pm.offset_cells[NPES_IO] * sizeof(cell_t);
  athread_join();
  timer_stop("import_box_sw");
  //printf("%ld %ld %ld %ld\n", ret, pm.offset_atoms[NPES_IO], pm.offset_cells[NPES_IO], estimate_atoms_size(fields));
  return 0;
}
#endif
