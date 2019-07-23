#include "data.h"
void esmd_pair_setup(esmd_t *md, ireal cutoff){
  md->pot_conf->cutoff = cutoff;
}

