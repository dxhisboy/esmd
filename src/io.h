#include <data.h>
enum rawio_fields {
  RAW_X = 1,
  RAW_V = 2,
  RAW_Q = 4,
  RAW_M = 8
};
void load_raw_x_atoms(esmd_t *md, const char *);
