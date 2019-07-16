#ifndef IO_H_
#define IO_H_
#include <data.h>
//function signatures
void load_raw_x_atoms(esmd_t *md, const char *path);
void load_raw_xv_atoms(esmd_t *md, const char *path);
void print_atoms_f(esmd_t *md);
//end function signatures
#endif
