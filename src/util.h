#ifndef UTIL_H_
#define UTIL_H_
#include <string.h>
#include <data.h>
enum transfer_flags {
  TRANS_ADJ_X = 1,
  TRANS_INC_F = 2,
  TRANS_TRUNC = 4
};
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define esmd_memcpy memcpy

//function signatures
size_t esmd_fields_size(int fields);
void esmd_export_cell(esmd_t *md, void *target, int fields, int flags, int celloff);
void esmd_import_cell(esmd_t *md, void *source, int fields, int flags, int celloff, areal *off);
size_t esmd_export_box(esmd_t *md, void *target, int fields, int flags, int xlo, int ylo, int zlo, int xlen, int ylen, int zlen);
size_t esmd_import_box(esmd_t *md, void *target, int fields, int flags,  int xlo, int ylo, int zlo, int xlen, int ylen, int zlen, areal *off);
//end function signatures
#endif
