#include <math.h>
#include "data.h"
#include "box.h"
#include "multiproc.h"
#include "pair.h"
#include "io.h"
#include "pair_lj.h"
#include "integrate.h"
#include <lattice.h>

int main(int argc, char **argv){
  esmd_t md;
  memory_init();
  MPI_Init(&argc, &argv);
  md.mpp.comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //esmd_pair_setup(&md, 2.5);
  areal cutoff = 2.5;
  ireal epsilon = 1.0;
  ireal sigma = 1.0;
  ireal mass = 1.0;
  md.integrate_conf.dt = 0.005;
  pair_lj_setup(&md, &cutoff, &epsilon, &sigma, &mass, 1);
  areal lg = 32 * pow((4.0 / 0.8442), (1.0 / 3.0));
  lattice_conf_t *lat_conf = &(md.lat_conf);
  lat_conf->atom_types = NULL;
  lat_conf->type = LAT_FCC;
  lat_conf->dens = 0.8442;
  lat_conf->nx = 32;
  lat_conf->ny = 32;
  lat_conf->nz = 32;
  
  md.utype = UNIT_LJ;

  esmd_lattice_scale(&md, lat_conf);
  esmd_set_box_size_by_lattice(&md);
  //esmd_set_box_size(&md, lg, lg, lg);
  esmd_box_setup_global(&md);
  esmd_multiproc_part_cart(&md, 2, 2, 1, rank);
  esmd_box_setup_local(&md);
  esmd_create_atoms_by_lattice(&md);
  thermo_init(&md);
  //printf("%d\n", md.natoms);
  scale_to_temp(&md, 1.44);

  esmd_exchange_cell(&md, LOCAL_TO_HALO, CELL_META | CELL_X | CELL_T, TRANS_ADJ_X);
  pair_lj_force(&md);
  esmd_exchange_cell(&md, HALO_TO_LOCAL, CELL_F, TRANS_INC_F);
  //return 0;
  for (int i = 0; i < 10; i ++){
    integrate(&md);
  }

  memory_print();
  MPI_Finalize();
  return 0;
}
