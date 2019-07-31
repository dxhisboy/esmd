#ifndef SWDATA_H_
#define SWDATA_H_
#include <data.h>
#define NCPE 64
#define BLK 4

#define __SW_CBC__(x) (x) * (x) * (x)
#define NNEIGHBOR ((__SW_CBC__(NCELL_CUT * 2 + 1) - 1) / 2 + 1)
typedef struct {
  //int start[64][3], count[64][3];
  areal (*frep)[CELL_SIZE][3];
  int read_from[NNEIGHBOR], store_to[NNEIGHBOR];
  int nupdate;
} swdata_t;
#endif
