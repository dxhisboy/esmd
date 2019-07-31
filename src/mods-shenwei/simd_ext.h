#include <simd.h>
#define vshuffd_rc(a, b, c, d) (d | (c << 2) | (b << 4) | (a << 6))

inline doublev4 simd_vsumd(doublev4 in){
  doublev4 ret = in;
  ret += simd_vshff(ret, ret, vshuffd_rc(2, 3, 0, 1));
  ret += simd_vshff(ret, ret, vshuffd_rc(1, 0, 3, 2));
  /* float tmp; */
  /* doublev4 ret; */
  /* //1,0,3,2 (01 00 11 10) 1011 0001 0xd1 */
  /* asm("vshff %2, %2, 0xb1, %1\n\t" */
  /*     "vaddd %2, %1, %1\n\t" */
  /*     "vshff %1, %1, 0x4e, %0\n\t" */
  /*     "vaddd %1, %0, %0\n\t" */
  /*     : "=r"(ret), "=r"(tmp) : "r"(in)); */
  return ret;
}
#define simd_vsumd_m(ret) {                                     \
    ret += simd_vshff(ret, ret, vshuffd_rc(2, 3, 0, 1));        \
    ret += simd_vshff(ret, ret, vshuffd_rc(1, 0, 3, 2));        \
  }
inline int simd_vmatchd(doublev4 in, int val){
  int ret;
  asm("vmatch %2, %1, %0\n\t"
      : "=r"(ret) : "r"(val), "r"(in));
  return ret;
}

#define simd_vmatchd_m(ret, in, val)   asm("vmatch %2, %1, %0\n\t" : "=r"(ret) : "r"(val), "r"(in));
#define simd_vselne(a, b, c) simd_vseleq(a, c, b)
inline void simd_load_4x3d(doublev4 *v0, doublev4 *v1, doublev4 *v2, void *st){
  doublev4 l0, l1, l2;
  simd_load(l0, st +  0);
  simd_load(l1, st + 32);
  simd_load(l2, st + 64);
  doublev4 t0 = simd_vshff(l2, l1, 0x6b);
  doublev4 t1 = simd_vshff(l1, l0, 0x46);
  *v0 = simd_vshff(t0, l0, 0xdc);
  *v1 = simd_vshff(t0, t1, 0x89);
  *v2 = simd_vshff(l2, t1, 0xcc);
}
inline void simd_cpyo(void *dest, void *src, size_t len){
  for (int i = 0; i < len; i += 128){
    uint256 tmp0, tmp1, tmp2, tmp3;
    simd_load(tmp0, src + i +  0);
    simd_load(tmp1, src + i + 32);
    simd_load(tmp2, src + i + 64);
    simd_load(tmp3, src + i + 96);
    
    simd_store(tmp0, dest + i +  0);
    simd_store(tmp1, dest + i + 32);
    simd_store(tmp2, dest + i + 64);
    simd_store(tmp3, dest + i + 96);
  }
}

inline void simd_zfillo(void *dest, size_t len){
  for (int i = 0; i < len; i += 128){
    asm("vstd $31,  0(%0)\n\t"
        "vstd $31, 32(%0)\n\t"
        "vstd $31, 64(%0)\n\t"
        "vstd $31, 96(%0)\n\t"
        : : "r"(dest + i));

  }
}
inline void simd_3x4_back(doublev4 *v0, doublev4 *v1, doublev4 *v2, doublev4 a0, doublev4 a1, doublev4 a2){
  doublev4 t0 = simd_vshff(a2, a1, 0x44);
  doublev4 t1 = simd_vshff(a2, a1, 0xee);
  doublev4 t2 = simd_vshff(a1, a0, 0xee);
  doublev4 t3 = simd_vshff(t0, a0, 0x84);
  doublev4 t4 = simd_vshff(t2, a2, 0xde);
  *v0 = simd_vshff(t3, t3, 0x78);
  *v1 = simd_vshff(t2, t0, 0x8d);
  *v2 = simd_vshff(t1, t4, 0xd8);
}

#define simd_3x4_back_m(v0, v1, v2, a0, a1, a2) {               \
    doublev4 t0, t1;                                            \
    asm("vshff %7, %6, 0x44, %0\n\t"                            \
        "vshff %7, %6, 0xee, %1\n\t"                            \
        "vshff %6, %5, 0xee, %2\n\t"                            \
        "vshff %0, %5, 0x84, %3\n\t"                            \
        "vshff %2, %7, 0xde, %4\n\t"                            \
        "vshff %3, %3, 0x78, %3\n\t"                            \
        "vshff %2, %0, 0x8d, %2\n\t"                            \
        "vshff %1, %4, 0xd8, %4\n\t"                            \
        : "=&r"(t0), "=&r"(t1), "=&r"(v1), "=&r"(v0), "=&r"(v2) \
        : "r"(a0), "r"(a1), "r"(a2));                           \
  }
//v0->t3, v1->t2, v2->t4
    /* doublev4 t0 = simd_vshff(a2, a1, 0x44);             \ */
    /* doublev4 t1 = simd_vshff(a2, a1, 0xee);             \ */
    /* doublev4 t2 = simd_vshff(a1, a0, 0xee);             \ */
    /* doublev4 t3 = simd_vshff(t0, a0, 0x84);             \ */
    /* doublev4 t4 = simd_vshff(t2, a2, 0xde);             \ */
    /* v0 = simd_vshff(t3, t3, 0x78);                     \ */
    /* v1 = simd_vshff(t2, t0, 0x8d);                     \ */
    /* v2 = simd_vshff(t1, t4, 0xd8);                     \ */
