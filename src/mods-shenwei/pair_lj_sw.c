/* #include <assert.h> */
/* #include <math.h> */
#include <geometry.h>
/* #include <timer.h> */

#define TEMPLATE <pair_lj_template_sw.h>
#ifdef MPE
#define FUNCTION slave_pair_lj_force_cpe
#endif
#include <ev_gen.h>
#ifdef MPE
#include <athread.h>
#include <timer.h>
void pair_lj_force_sw(esmd_t *md, int evflag){
  timer_start("force_sw");
  athread_spawn(pair_lj_force_cpe_vers[evflag], md);
  athread_join();
  timer_stop("force_sw");
}
#endif
