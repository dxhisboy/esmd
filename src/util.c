#include <data.h>
#include <util.h>
#include <timer.h>
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
#define export_field(fcode, fname, size)			\
  if (fields & fcode) {						\
    esmd_memcpy(buffer + offset, data->fname + start, (size));	\
    offset += size;						\
  }

inline int export_data(void *buffer, cell_t *cell, celldata_t *data,
			int start, int end, int fields, int flags){
  int offset = 0;
  int natoms = end - start;
  export_field(CELL_X, x, sizeof(areal) * natoms * 3);
  export_field(CELL_V, v, sizeof(areal) * natoms * 3);
  export_field(CELL_F, f, sizeof(areal) * natoms * 3);
  export_field(CELL_Q, q, sizeof(areal) * natoms);
  export_field(CELL_T, type, sizeof(int) * natoms);
  export_field(CELL_E, export, sizeof(int) * natoms);
  return offset;
}
int esmd_export_cell(esmd_t *md, void *buffer, int fields, int flags, int celloff){
  cell_t *cell = md->box.cells + celloff;
  celldata_t *data = md->box.celldata + celloff;

  size_t offset = 0;
  if (fields & CELL_META) {
    esmd_memcpy(buffer, md->box.cells + celloff, sizeof(cell_t));
    offset += sizeof(cell_t);
  }

  if (flags & (TRANS_ATOMS | TRANS_EXPORTS)) {
    if (flags & TRANS_ATOMS){
      offset += export_data(buffer + offset, cell, data, 0, cell->natoms, fields, flags);
    }
    if (flags & TRANS_EXPORTS){
      offset += export_data(buffer + offset, cell, data, cell->export_ptr, CELL_SIZE, fields, flags);
    }
  } else {
    offset += export_data(buffer + offset, cell, data, 0, CELL_SIZE, fields, flags);
  }
  return offset;
}

#define import_field(fcode, fname, size)			\
  if (fields & fcode) {						\
    esmd_memcpy(data->fname + start, buffer + offset, (size));	\
    offset += size;						\
  }

inline int import_data(void *buffer, cell_t *cell, celldata_t *data,
			int start, int end, int fields, int flags, areal *off){
  size_t offset = 0;
  int natoms = end - start;
  import_field(CELL_X, x, sizeof(areal) * natoms * 3);
  import_field(CELL_V, v, sizeof(areal) * natoms * 3);
  if (CELL_F & fields){
    if (TRANS_INC_F & flags) {
      areal (*cellf)[3] = data->f;
      areal (*importf)[3] = buffer + offset;
      for (int i = 0; i < cell->natoms; i ++){
        cellf[start + i][0] += importf[i][0];
        cellf[start + i][1] += importf[i][1];
        cellf[start + i][2] += importf[i][2];
      }
      offset += sizeof(areal) * natoms * 3;
    } else {
      import_field(CELL_F, f, sizeof(areal) * natoms * 3);
    }
  }

  import_field(CELL_Q, q, sizeof(areal) * natoms);
  import_field(CELL_T, type, sizeof(int) * natoms);
  import_field(CELL_E, export, sizeof(int) * natoms);
  if ((TRANS_ADJ_X & flags) && (CELL_X & fields)) {
    areal (*cellx)[3] = data->x;
    for (int i = 0; i < natoms; i ++){
      cellx[start + i][0] += off[0];
      cellx[start + i][1] += off[1];
      cellx[start + i][2] += off[2];
    }
  }
  return offset;
}
int esmd_import_cell(esmd_t *md, void *buffer, int fields, int flags, int celloff, areal *off){
  cell_t *cell = md->box.cells + celloff;
  celldata_t *data = md->box.celldata + celloff;
  size_t offset = 0;
  if (fields & CELL_META) {
    esmd_memcpy(md->box.cells + celloff, buffer, sizeof(cell_t));
    offset += sizeof(cell_t);
    md->box.cells[celloff].bbox_ideal[0][0] += off[0];
    md->box.cells[celloff].bbox_ideal[0][1] += off[1];
    md->box.cells[celloff].bbox_ideal[0][2] += off[2];
    md->box.cells[celloff].bbox_ideal[1][0] += off[0];
    md->box.cells[celloff].bbox_ideal[1][1] += off[1];
    md->box.cells[celloff].bbox_ideal[1][2] += off[2];
  }
  if (flags & (TRANS_ATOMS | TRANS_EXPORTS)) {
    if (flags & TRANS_ATOMS){
      offset += import_data(buffer + offset, cell, data, 0, cell->natoms, fields, flags, off);
    }
    if (flags & TRANS_EXPORTS){
      offset += import_data(buffer + offset, cell, data, cell->export_ptr, CELL_SIZE, fields, flags, off);
    }
  } else {
    offset += import_data(buffer + offset, cell, data, 0, CELL_SIZE, fields, flags, off);
  }
  return offset;
}

size_t esmd_export_box(esmd_t *md, void *buffer, int fields, int flags, int xlo, int ylo, int zlo, int xlen, int ylen, int zlen){
  timer_start("esmd_export_box");
  box_t *box = &(md->box);
  size_t entry_size = esmd_fields_size(fields);
  size_t bufoff = 0;
  for (int k = zlo; k < zlo + zlen; k ++){
    for (int j = ylo; j < ylo + ylen; j ++){
      for (int i = xlo; i < xlo + xlen; i ++){
        int celloff = get_cell_off(box, i, j, k);
        bufoff += esmd_export_cell(md, buffer + bufoff, fields, flags, celloff);
      }
    }
  }
  timer_stop("esmd_export_box");
  return bufoff;
  //return entry_size * xlen * ylen * zlen;
}

size_t esmd_import_box(esmd_t *md, void *buffer, int fields, int flags,  int xlo, int ylo, int zlo, int xlen, int ylen, int zlen, areal *off){
  timer_start("esmd_import_box");
  box_t *box = &(md->box);
  int entry_size = esmd_fields_size(fields);
  size_t bufoff = 0;
  for (int k = zlo; k < zlo + zlen; k ++){
    for (int j = ylo; j < ylo + ylen; j ++){
      for (int i = xlo; i < xlo + xlen; i ++){
        int celloff = get_cell_off(box, i, j, k);
        bufoff += esmd_import_cell(md, buffer + bufoff, fields, flags, celloff, off);
      }
    }
  }
  timer_stop("esmd_import_box");
  return bufoff;
  //return entry_size * xlen * ylen * zlen;
}
