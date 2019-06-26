#include <data.h>
#include <util.h>
#define add_size(fcode, size) if (fields & fcode) accu_size += (size);
size_t esmd_fields_size(int fields) {
  size_t accu_size = 0;
  add_size(CELL_META, sizeof(cell_t));
  add_size(CELL_X, sizeof(areal) * CELL_SIZE * 3);
  add_size(CELL_V, sizeof(areal) * CELL_SIZE * 3);
  add_size(CELL_F, sizeof(areal) * CELL_SIZE * 3);
  add_size(CELL_Q, sizeof(areal) * CELL_SIZE);
  add_size(CELL_T, sizeof(int) * CELL_SIZE);
  add_size(CELL_E, sizeof(int) * CELL_SIZE);
  return accu_size;
}
#define export_field(fcode, fname, size)                                \
  if (fields & fcode) {                                                 \
    esmd_memcpy(target + nbytes_write, md->box.celldata[celloff].fname, (size)); \
    nbytes_write += size;                                               \
  }

void esmd_export_cell(esmd_t *md, void *target, int fields, int flags, int celloff){
  size_t nbytes_write = 0;
  if (fields & CELL_META) {
    esmd_memcpy(target, md->box.cells + celloff, sizeof(cell_t));
    nbytes_write += sizeof(cell_t);
  }
  export_field(CELL_X, x, sizeof(areal) * CELL_SIZE * 3);
  export_field(CELL_V, v, sizeof(areal) * CELL_SIZE * 3);
  export_field(CELL_Q, q, sizeof(areal) * CELL_SIZE);
  export_field(CELL_T, type, sizeof(int) * CELL_SIZE);
  export_field(CELL_F, f, sizeof(areal) * CELL_SIZE * 3);
  export_field(CELL_E, export, sizeof(int) * CELL_SIZE);
  if ((CELL_F & fields) && (TRANS_INC_F & flags)) {
    memset(md->box.celldata[celloff].f, 0, sizeof(areal) * CELL_SIZE * 3);
  }
}

#define import_field(fcode, fname, size)                                \
  if (fields & fcode) {                                                 \
    esmd_memcpy(md->box.celldata[celloff].fname, source + nbytes_read, (size)); \
    nbytes_read += size;                                                \
  }

void esmd_import_cell(esmd_t *md, void *source, int fields, int flags, int celloff, areal *off){
  size_t nbytes_read = 0;
  if (fields & CELL_META) {
    esmd_memcpy(md->box.cells + celloff, source, sizeof(cell_t));
    nbytes_read += sizeof(cell_t);
    md->box.cells[celloff].bbox_ideal[0][0] += off[0];
    md->box.cells[celloff].bbox_ideal[0][1] += off[1];
    md->box.cells[celloff].bbox_ideal[0][2] += off[2];
    md->box.cells[celloff].bbox_ideal[1][0] += off[0];
    md->box.cells[celloff].bbox_ideal[1][1] += off[1];
    md->box.cells[celloff].bbox_ideal[1][2] += off[2];
  }
  import_field(CELL_X, x, sizeof(areal) * CELL_SIZE * 3);
  import_field(CELL_V, v, sizeof(areal) * CELL_SIZE * 3);
  import_field(CELL_Q, q, sizeof(areal) * CELL_SIZE);
  import_field(CELL_T, type, sizeof(int) * CELL_SIZE);
  if (CELL_F & fields){
    if (TRANS_INC_F & flags) {
      areal (*cellf)[3] = md->box.celldata[celloff].f;
      areal (*importf)[3] = source + nbytes_read;
      for (int i = 0; i < md->box.cells[celloff].natoms; i ++){
        cellf[i][0] += importf[i][0];
        cellf[i][1] += importf[i][1];
        cellf[i][2] += importf[i][2];
      }
      nbytes_read += sizeof(areal) * CELL_SIZE * 3;
    } else {
      import_field(CELL_F, f, sizeof(areal) * CELL_SIZE * 3);
    }
  }
  import_field(CELL_E, export, sizeof(int) * CELL_SIZE);
  if ((TRANS_ADJ_X & flags) && (CELL_X & fields)) {
    areal (*cellx)[3] = md->box.celldata[celloff].x;
    for (int i = 0; i < CELL_SIZE; i ++){
      cellx[i][0] += off[0];
      cellx[i][1] += off[1];
      cellx[i][2] += off[2];
    }
  }
}

size_t esmd_export_box(esmd_t *md, void *target, int fields, int flags, int xlo, int ylo, int zlo, int xlen, int ylen, int zlen){
  box_t *box = &(md->box);
  size_t entry_size = esmd_fields_size(fields);

  for (int i = xlo; i < xlo + xlen; i ++){
    for (int j = ylo; j < ylo + ylen; j ++){
      for (int k = zlo; k < zlo + zlen; k ++){
        int bufoff = (((i - xlo) * ylen + j - ylo) * zlen + k - zlo) * entry_size;
        int celloff = get_cell_off(box, i, j, k);
        esmd_export_cell(md, target + bufoff, fields, flags, celloff);
      }
    }
  }
  return entry_size * xlen * ylen * zlen;
}

size_t esmd_import_box(esmd_t *md, void *target, int fields, int flags,  int xlo, int ylo, int zlo, int xlen, int ylen, int zlen, areal *off){
  box_t *box = &(md->box);
  int entry_size = esmd_fields_size(fields);

  for (int i = xlo; i < xlo + xlen; i ++){
    for (int j = ylo; j < ylo + ylen; j ++){
      for (int k = zlo; k < zlo + zlen; k ++){
        int bufoff = (((i - xlo) * ylen + j - ylo) * zlen + k - zlo) * entry_size;
        int celloff = get_cell_off(box, i, j, k);
        esmd_import_cell(md, target + bufoff, fields, flags, celloff, off);
      }
    }
  }
  return entry_size * xlen * ylen * zlen;
}
