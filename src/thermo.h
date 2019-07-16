#ifndef THERMO_H_
#define THERMO_H_
#include <data.h>
//function signatures
double temperature(esmd_t *md);
double compute_kinetic_local(esmd_t *md);
void thermo_init(esmd_t *md);
void thermo_compute(esmd_t *md);
void scale_to_temp(esmd_t *md, areal t_req);
//end function signatures
#endif
