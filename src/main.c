#include <math.h>
#include "data.h"
#include "box.h"
#include "multiproc.h"
#include "pair.h"
#include "io.h"
#include "pair_lj.h"
#include "integrate.h"
#include <timer.h>
#include <thermo.h>
#include <lattice.h>
#include <log.h>
int main(int argc, char **argv){
  esmd_t md;
  memory_init();
  timer_init();
  esmd_mpi_init(&md);
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
  lat_conf->nx = 128;
  lat_conf->ny = 32;
  lat_conf->nz = 32;
  
  md.utype = UNIT_LJ;

  esmd_lattice_scale(&md, lat_conf);
  esmd_set_box_size_by_lattice(&md);
  esmd_box_setup_global(&md);
  esmd_auto_part(&md);
  esmd_multiproc_part(&md);
  master_info("divided tasks to %d*%d*%d grid\n", md.mpp.npx, md.mpp.npy, md.mpp.npz);
  esmd_box_setup_local(&md);
  esmd_create_atoms_by_lattice(&md);
  report_cell_info(&md);
  thermo_init(&md);
  scale_to_temp(&md, 1.44);

  esmd_exchange_cell(&md, LOCAL_TO_HALO, CELL_META | CELL_X | CELL_T, TRANS_ADJ_X);
  pair_lj_force(&md);
  esmd_exchange_cell(&md, HALO_TO_LOCAL, CELL_F, TRANS_INC_F);
  //return 0;
  for (int i = 0; i < 10; i ++){
    timer_start("integrate");
    integrate(&md);
    timer_stop("integrate");
  }
  report_cell_info(&md);

  //memory_print();
  timer_print(md.mpp.comm);
  MPI_Finalize();
  return 0;
}
