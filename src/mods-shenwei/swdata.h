#ifndef SWDATA_H_
#define SWDATA_H_
#include <data.h>
#define NCPE 64
#define BLK 4

typedef struct {
  int start[64][3], count[64][3];
  areal (*frep)[CELL_SIZE][3];
  int nupdate;
} swdata_t;
#endif
