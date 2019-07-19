#define __CAT__(x, y) x ## _ ## y
#define CAT(x, y) __CAT__(x, y)

#define VIRIAL_CODE 1
#define ENERGY_CODE 2
#define VER_CODE 0
#include <ev_ver.h>
#undef VER_CODE

#define VER_CODE 1
#include <ev_ver.h>
#undef VER_CODE

#define VER_CODE 2
#include <ev_ver.h>
#undef VER_CODE

#define VER_CODE 3
#include <ev_ver.h>
#undef VER_CODE

#ifdef FUNCTION
void (*pair_lj_force_vers[4])(esmd_t *) = {
  CAT(FUNCTION, 0),
  CAT(FUNCTION, 1),
  CAT(FUNCTION, 2),
  CAT(FUNCTION, 3)
};
#endif
#undef __CAT__
#undef CAT
