#include <data.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <geometry.h>
#include <timer.h>
//#define DEBUG_THIS_FILE
#include <log.h>

void pair_lj_setup(esmd_t *md, areal *cutoff, ireal *epsilon, ireal *sigma, ireal *mass, int ntypes){
  assert(ntypes < MAX_TYPES);
  md->pot_conf->cutoff = 0;
  md->pot_conf->ntypes = ntypes;
  for (int i = 0; i < ntypes; i ++){
    md->pot_conf->rmass[i] = 1.0 / mass[i];
    md->pot_conf->mass[i] = mass[i];
    for (int j = 0; j < ntypes; j ++){
      ireal sigma2 = sigma[i * ntypes + j] * sigma[i * ntypes + j];
      ireal sigma6 = sigma2 * sigma2 * sigma2;
      ireal sigma12 = sigma6 * sigma6;
      md->pot_conf->param.lj.c6[i][j] = 4.0 * 6.0 * epsilon[i * ntypes + j] * sigma6;
      md->pot_conf->param.lj.c12[i][j] = 4.0 * 12.0 * epsilon[i * ntypes + j] * sigma12;
      md->pot_conf->param.lj.ec6[i][j] = 4.0 * epsilon[i * ntypes + j] * sigma6;
      md->pot_conf->param.lj.ec12[i][j] = 4.0 * epsilon[i * ntypes + j] * sigma12;
      md->pot_conf->param.lj.cutoff2[i][j] = cutoff[i * ntypes + j] * cutoff[i * ntypes + j];
      if (cutoff[i * ntypes + j] > md->pot_conf->cutoff){
        md->pot_conf->cutoff = cutoff[i * ntypes + j];
      }
    }
  }
}

#define TEMPLATE <pair_lj_template.h>
#define FUNCTION pair_lj_force
#include <ev_gen.h>
void pair_lj_force(esmd_t *md, int evflag) {
  timer_start("force");
  pair_lj_force_vers[evflag](md);
  timer_stop("force");
#ifdef DEBUG_THIS_FILE
  areal evdwl_gbl;
  esmd_global_sum_scalar(md, &evdwl_gbl, evdwl);

  if (md->mpp->pid == 0){
    debug("%d %f %d %d\n", md->mpp->pid, evdwl_gbl, total_atoms, total_int);
  }
#endif
}
