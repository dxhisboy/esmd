/* #include <assert.h> */
/* #include <math.h> */
#include <geometry.h>
/* #include <timer.h> */

#ifdef CPE
#include <simd.h>
#include <slave.h>
#include <dma_macros.h>
#include <reg_reduce.h>
#include <swdata.h>
#include <simd_ext.h>
#define LWPF_UNIT U(PAIR_LJ)
#define LWPF_KERNELS K(ALL) K(FILL) K(COMP) K(SYN) K(FINI) K(RI) K(RJ) K(WJ) K(B2B) K(CONV) K(PROF)
#undef inline
#include <lwpf2/lwpf2.h>
#endif

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
