#include <data.h>
void camd_pair_setup(camd_t *camd, ireal cutoff){
  camd->pair_conf.cutoff = cutoff;
}

