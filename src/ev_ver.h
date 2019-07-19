#ifdef VIRIAL
#undef VIRIAL
#endif

#if (VER_CODE & VIRIAL_CODE)
#define VIRIAL 1
#else
#define VIRIAL 0
#endif


#ifdef ENERGY
#undef ENERGY
#endif

#if (VER_CODE & ENERGY_CODE)
#define ENERGY 1
#else
#define ENERGY 0
#endif

#include TEMPLATE
