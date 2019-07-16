#ifndef BOX_H_
#define BOX_H_
#include <data.h>
//function signatures
void esmd_set_box_size(esmd_t *md, areal x, areal y, areal z);
void esmd_box_setup_global(esmd_t *md);
void esmd_box_setup_local(esmd_t *md);
void box_add_atom(box_t *box, areal *x, areal *v, ireal q, int type);
//end function signatures
#endif
