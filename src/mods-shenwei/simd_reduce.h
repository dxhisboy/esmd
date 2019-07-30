#include <simd.h>
inline doublev4 simd_vsumd(doublev4 in){
  float tmp;
  doublev4 ret;
  //1,0,3,2 (01 00 11 10) 1011 0001 0xd1
  asm("vshff %2, %2, 0xb1, %1\n\t"
      "vaddd %2, %1, %1\n\t"
      "vshff %1, %1, 0x4e, %0\n\t"
      "vaddd %1, %0, %0\n\t"
      : "=r"(ret), "=r"(tmp) : "r"(in));
  return ret;
}
