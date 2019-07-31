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
#define DEBUG_THIS_FILE
#include <log.h>
#include <swlu.h>
int main(int argc, char **argv){
  esmd_t md;
  memory_init();
  md.box = esmd_malloc(sizeof(box_t), "box meta");
  md.mpp = esmd_malloc(sizeof(multiproc_t), "mpi meta");
  md.pot_conf = esmd_malloc(sizeof(potential_conf_t), "potential config");

  timer_init();
  esmd_mpi_init(&md);
  //esmd_pair_setup(&md, 2.5);
  areal cutoff = 2.5;
  ireal epsilon = 1.0;
  ireal sigma = 1.0;
  ireal mass = 1.0;
  md.dt = 0.005;
  pair_lj_setup(&md, &cutoff, &epsilon, &sigma, &mass, 1);
  areal lg = 32 * pow((4.0 / 0.8442), (1.0 / 3.0));
  lattice_conf_t *lat_conf = &(md.lat_conf);
  lat_conf->atom_types = NULL;
  lat_conf->type = LAT_FCC;
  lat_conf->dens = 0.8442;
  lat_conf->nx = 80;
  lat_conf->ny = 40;
  lat_conf->nz = 40;
  
  md.utype = UNIT_LJ;

  esmd_lattice_scale(&md, lat_conf);

  debug("sizeof md: %ld\n", sizeof(md));
  debug("sizeof box: %ld\n", sizeof(md.box));
  debug("sizeof pot_conf: %ld\n", sizeof(md.pot_conf));
  esmd_set_box_size_by_lattice(&md);
  master_info("box size is %f %f %f\n", md.box->lglobal[0], md.box->lglobal[1], md.box->lglobal[2]);

  master_info("lattice scale is %f\n", lat_conf->scale);
  esmd_box_setup_global(&md);
  master_info("cell size is %f %f %f\n", md.box->lcell[0], md.box->lcell[1], md.box->lcell[2]);
  esmd_auto_part(&md);
  esmd_multiproc_part(&md);
  master_info("divided tasks to %d*%d*%d grid\n", md.mpp->npx, md.mpp->npy, md.mpp->npz);
  esmd_box_setup_local(&md);
  esmd_create_atoms_by_lattice(&md);
  report_cell_info(&md);
  thermo_init(&md);
  scale_to_temp(&md, 1.44);
  md.accu_local.virial = 0;
  md.accu_local.epot = 0;
  md.accu_local.kinetic = 0;
  md.nthermo = 100;
  shenwei_init(&md);
  esmd_exchange_cell(&md, LOCAL_TO_HALO, CELL_META | CELL_X | CELL_T, TRANS_ADJ_X | TRANS_ATOMS);
  pair_lj_force(&md, 3);
  esmd_exchange_cell(&md, HALO_TO_LOCAL, CELL_F, TRANS_INC_F | TRANS_ATOMS);
  //return 0;
  md.step = 1;
  /* swlu_prof_init(); */
  /* swlu_prof_start(); */
  //enable_memcpy_log();
  for (int i = 0; i < 1000; i ++){
    timer_start("integrate");
    integrate(&md);
    timer_stop("integrate");
  }
  //disable_memcpy_log();
  /* swlu_prof_stop(); */
  /* swlu_prof_print(); */
  report_cell_info(&md);

  timer_print(md.mpp->comm);
  shenwei_destroy();
  MPI_Finalize();
  return 0;
}
