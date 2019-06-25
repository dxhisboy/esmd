#include "data.h"
#include "box.h"
#include "multiproc.h"
#include "pair.h"
#include "io.h"
#include "pair_lj.h"
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
  pair_lj_setup(&md, &cutoff, &epsilon, &sigma, 1);
  esmd_box_setup_global(&md, 53.747078, 53.747078, 53.747078);
  esmd_multiproc_part_cart(&md, 1, 1, 1, rank);
  esmd_box_setup_local(&md);
  load_raw_x_atoms(&md, "/home/xduan/workspace/miniMD/ref/data/x.bin");
  esmd_exchange_cell(&md, CELL_META | CELL_X | CELL_T);
  pair_lj_force(&md);
  print_atoms(&md);
  MPI_Finalize();
  return 0;
}
