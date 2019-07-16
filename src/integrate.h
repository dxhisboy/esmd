#ifndef INTEGRATE_H_
#define INTEGRATE_H_
#include <data.h>
//function signatures
void initial_integrate_nve(esmd_t *md);
void final_integrate_nve(esmd_t *md);
void esmd_export_atoms(esmd_t *md);
int esmd_import_atoms_pairwise(box_t *box, int self_off, int neigh_off);
void esmd_import_atoms_outer_from_local(esmd_t *md);
void esmd_import_atoms_outer_from_halo(esmd_t *md);
void esmd_import_atoms_inner_from_local(esmd_t *md);
int get_halo_type(box_t *box, int ii, int jj, int kk);
void esmd_import_atoms_halo_from_local(esmd_t *md);
void esmd_import_atoms(esmd_t *md);
void integrate(esmd_t *md);
//end function signatures
#endif
