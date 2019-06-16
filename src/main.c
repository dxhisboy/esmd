#include "data.h"
#include "box.h"
#include "multiproc.h"
int main(int argc, char **argv){
  esmd_t md;
  memory_init();
  esmd_pair_setup(&md, 2.5);
  esmd_box_setup_global(&md, 32, 32, 32);
  esmd_multiproc_part_cart(&md, 1, 1, 1, 0);
  esmd_box_setup_local(&md);
  load_raw_x_atoms(&md, "/home/xduan/workspace/miniMD/ref/data/x.bin");
  
}
