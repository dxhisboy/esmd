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
    esmd_memcpy(target + nbytes_write, data->fname + start, (size)); \
    nbytes_write += size;                                               \
  }

inline size_t export_data(void *target, cell_t *cell, celldata_t *data,
			  int start, int end, int fields, int flags){
  int natoms = end - start;
  size_t nbytes_write = 0;

  export_field(CELL_X, x     , sizeof(areal) * natoms * 3);
  export_field(CELL_V, v     , sizeof(areal) * natoms * 3);
  export_field(CELL_F, f     , sizeof(areal) * natoms * 3);
  export_field(CELL_Q, q     , sizeof(ireal) * natoms    );
  export_field(CELL_T, type  , sizeof(int)   * natoms    );
  export_field(CELL_E, export, sizeof(int)   * natoms    );
}

void esmd_export_cell(esmd_t *md, void *target, int fields, int flags, int celloff){
  size_t nbytes_write = 0;
  cell_t *cell = md->box.cells + celloff;
  celldata_t *celldata = md->box.celldata + celloff;
  if (fields & CELL_META) {
    esmd_memcpy(target, cell, sizeof(cell_t));
    nbytes_write += sizeof(cell_t);
  }
  int natom_trans = CELL_SIZE;
  export_data(target + nbytes_write, cell, celldata, 0, natom_trans, fields, flags);
  /* if (flags & TRANS_ATOMONLY) { */
  /*   natom_trans = md->box.cells[celloff].natoms; */
  /* } */
  /* if (flags & TRANS_ATOMEXPORT) { */
  /*   natom_trans = md->box.cells[celloff].natoms + md->box.cells[celloff].natoms; */
  /* } */

  /* export_field(CELL_X, x, sizeof(areal) * natom_trans * 3); */
  /* export_field(CELL_V, v, sizeof(areal) * natom_trans * 3); */
  /* export_field(CELL_F, f, sizeof(areal) * natom_trans * 3); */
  /* export_field(CELL_Q, q, sizeof(ireal) * natom_trans); */
  /* export_field(CELL_T, type, sizeof(int) * natom_trans); */
  /* export_field(CELL_E, export, sizeof(int) * natom_trans); */
  /* if ((CELL_F & fields) && (TRANS_INC_F & flags)) { */
  /*   memset(md->box.celldata[celloff].f, 0, sizeof(areal) * natom_trans * 3); */
  /* } */
}

/* #define import_field(fcode, fname, size)                                \ */
/*   if (fields & fcode) {                                                 \ */
/*     esmd_memcpy(md->box.celldata[celloff].fname, source + nbytes_read, (size)); \ */
/*     nbytes_read += size;                                                \ */
/*   } */

#define import_field(fcode, fname, size)                                \
  if (fields & fcode) {                                                 \
    esmd_memcpy(data->fname + start, source + nbytes_read, (size));	\
    nbytes_read += size;                                                \
  }

inline size_t esmd_import_data(void *source, cell_t *cell, celldata_t *data,
			       int start, int end, int fields, int flags, areal *off){
  int natoms = end - start;
  size_t nbytes_read = 0;
  import_field(CELL_X, x     , sizeof(areal) * natoms * 3);
  import_field(CELL_V, v     , sizeof(areal) * natoms * 3);
  if ((CELL_F & fields) && (flags & TRANS_INC_F)) {
    areal (*f)[3] = data->f;
    areal (*f_import)[3] = source + nbytes_read;
    for (int i = 0; i < natoms; i ++){
      f[start + i][0] = f_import[i][0];
      f[start + i][1] = f_import[i][1];
      f[start + i][2] = f_import[i][2];		  
    }
  } else {
    import_field(CELL_F, f     , sizeof(areal) * natoms * 3);
  }
  import_field(CELL_Q, q     , sizeof(ireal) * natoms    );
  import_field(CELL_T, type  , sizeof(int)   * natoms    );
  import_field(CELL_E, export, sizeof(int)   * natoms    );
  
  if ((TRANS_ADJ_X & flags) && (CELL_X & fields)) {
    areal (*x)[3] = data->x;
    for (int i = 0; i < natoms; i ++){
      x[start + i][0] += off[0];
      x[start + i][1] += off[1];
      x[start + i][2] += off[2];
    }
  }

}
void esmd_import_cell(esmd_t *md, void *source, int fields, int flags, int celloff, areal *off){
  cell_t *cell = md->box.cells + celloff;
  celldata_t *data = md->box.celldata + celloff;
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
  int natom_trans = CELL_SIZE;
  esmd_import_data(source + nbytes_read, cell, data, 0, natom_trans, fields, flags, off);
}

size_t esmd_export_box(esmd_t *md, void *target, int fields, int flags, int xlo, int ylo, int zlo, int xlen, int ylen, int zlen){
  box_t *box = &(md->box);
  size_t entry_size = esmd_fields_size(fields);
  for (int k = zlo; k < zlo + zlen; k ++){
    for (int j = ylo; j < ylo + ylen; j ++){
      for (int i = xlo; i < xlo + xlen; i ++){
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
  for (int k = zlo; k < zlo + zlen; k ++){
    for (int j = ylo; j < ylo + ylen; j ++){
      for (int i = xlo; i < xlo + xlen; i ++){
        int bufoff = (((i - xlo) * ylen + j - ylo) * zlen + k - zlo) * entry_size;
        int celloff = get_cell_off(box, i, j, k);
        esmd_import_cell(md, target + bufoff, fields, flags, celloff, off);
      }
    }
  }
  return entry_size * xlen * ylen * zlen;
}
