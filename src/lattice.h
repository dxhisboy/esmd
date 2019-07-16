#ifndef LATTICE_H_
#define LATTICE_H_
#include <data.h>
//function signatures
double minimd_random(int *seed);
double next_rand(int *seed);
void esmd_lattice_scale(esmd_t *md, lattice_conf_t *conf);
void esmd_set_box_size_by_lattice(esmd_t *md);
void esmd_create_atoms_by_lattice(esmd_t *md);
//end function signatures
#endif
