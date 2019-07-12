#include "data.h"
#include "box.h"
#include "multiproc.h"
#include "pair.h"
#include "io.h"
#include "pair_lj.h"
#include <math.h>
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
  esmd_box_setup_global(&md, lg, lg, lg);
  esmd_multiproc_part_cart(&md, 1, 1, 1, rank);
  esmd_box_setup_local(&md);
  load_raw_xv_atoms(&md, "/uni-mainz.de/homes/xiaoduan/miniMD/ref/data/xv.bin");
  //esmd_exchange_cell_local_to_halo(&md, CELL_META | CELL_X | CELL_T, TRANS_ADJ_X);
  //esmd_exchange_cell(&md, LOCAL_TO_HALO, CELL_META | CELL_X | CELL_T, TRANS_ADJ_X);
  esmd_exchange_cell(&md, LOCAL_TO_HALO, CELL_META | CELL_X | CELL_T, TRANS_ADJ_X);
  pair_lj_force(&md);
  esmd_exchange_cell(&md, HALO_TO_LOCAL, CELL_F, TRANS_INC_F);
  for (int i = 0; i < 10; i ++){
    integrate(&md);
  }

  memory_print();
  MPI_Finalize();
  return 0;
}
