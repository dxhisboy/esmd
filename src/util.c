#include <data.h>
#define add_size(fcode, size) if (fields & fcode) accu_size += (size);
int esmd_fields_size(int fields) {
  size_t accu_size = 0;
  add_size(CELL_META, sizeof(cell_t));
  add_size(CELL_X, sizeof(areal) * CELL_SIZE * 3);
  add_size(CELL_V, sizeof(areal) * CELL_SIZE * 3);
  add_size(CELL_F, sizeof(areal) * CELL_SIZE * 3);
  add_size(CELL_Q, sizeof(areal) * CELL_SIZE);
  add_size(CELL_T, sizeof(int) * CELL_SIZE);
  return accu_size;
}
#define export_field(fcode, fname, size)                                \
  if (fields & fcode) {                                                 \
    esmd_memcpy(target + nbytes_write, md->box->celldata[celloff].fname, (size)); \
    nbytes_write += size;                                               \
  }

void esmd_export_cell(esmd_t *md, void *target, int fields, int celloff){
  size_t nbytes_write = 0;
  if (fields & CELL_META) {
    esmd_memcpy(target, md->box->cell + celloff, sizeof(cell_t));
    nbytes_write += sizeof(cell_t);
  }
  export_field(CELL_X, x, sizeof(areal) * CELL_SIZE * 3);
  export_field(CELL_V, v, sizeof(areal) * CELL_SIZE * 3);
  export_field(CELL_F, f, sizeof(areal) * CELL_SIZE * 3);
  export_field(CELL_Q, q, sizeof(areal) * CELL_SIZE);
  export_field(CELL_T, type, sizeof(int) * CELL_SIZE);
}

#define import_field(fcode, fname, size)                                \
  if (fields & fcode) {                                                 \
    esmd_memcpy(md->box->celldata[celloff].fname, source + nbytes_write, (size)); \
    nbytes_write += size;                      \
  }

void esmd_import_cell(esmd_t *md, void *source, int fields, int celloff){
  size_t nbytes_read = 0;
  if (fields & CELL_META) {
    esmd_memcpy(md->box->cell + celloff, target, sizeof(cell_t));
    nbytes_read += sizeof(cell_t);
  }
  import_field(CELL_X, x, sizeof(areal) * CELL_SIZE * 3);
  import_field(CELL_V, v, sizeof(areal) * CELL_SIZE * 3);
  import_field(CELL_F, f, sizeof(areal) * CELL_SIZE * 3);
  import_field(CELL_Q, q, sizeof(areal) * CELL_SIZE);
  import_field(CELL_T, type, sizeof(int) * CELL_SIZE);
}
