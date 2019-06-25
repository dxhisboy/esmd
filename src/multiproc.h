#ifndef MULTIPROC_H_
#define MULTIPROC_H_
#include <util.h>
enum exchange_direction {
  LOCAL_TO_HALO = 0,
  HALO_TO_LOCAL = 1
};
void esmd_multiproc_part_cart(esmd_t *, int, int, int, int);
//void esmd_exchange_cell(esmd_t *md, int direction, int fields, int flags);
#endif
