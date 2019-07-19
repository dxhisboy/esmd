#ifndef PAIR_LJ_H_
#define PAIR_LJ_H_
#include <data.h>
//function signatures
void pair_lj_setup(esmd_t *md, areal *cutoff, ireal *epsilon, ireal *sigma, ireal *mass, int ntypes);
void pair_lj_force(esmd_t *md, int evflag);
//end function signatures
#endif
