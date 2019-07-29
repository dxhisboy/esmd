#include <stdio.h>
#include <string.h>
#include <simd.h>

typedef struct {
  void *dest;
  void *src;
  size_t cnt;
} memcpy_arg_t;
#ifdef CPE
#define BLKSZ 0x400L
#include <dma_macros.h>
void memcpy_cpe(memcpy_arg_t *gl_arg){
  if (_COL) return;
  dma_init();
  memcpy_arg_t arg;
  pe_get(gl_arg, &arg, sizeof(memcpy_arg_t));
  dma_syn();
  void *dest = arg.dest;
  void *src = arg.src;
  size_t cnt = arg.cnt;
  char buffer[BLKSZ];
  for (size_t blkst = BLKSZ * _ROW; blkst < cnt; blkst += BLKSZ * 8){
    int blksz = BLKSZ;
    if (blkst + blksz > cnt) blksz = cnt - blkst;
    pe_get(src + blkst, buffer, blksz);
    dma_syn();
    pe_put(dest + blkst, buffer, blksz);
    dma_syn();
  }
  
}
#endif
#ifdef MPE
inline int check_align(long val){
  return (val & -val);
}
void *swmemcpy_vec4(void *dest, void *src, size_t cnt){
  for (int i = 0; i < cnt; i += 128){
    intv8 tmp0, tmp1, tmp2, tmp3;
    asm("vldd %0,  0(%5)\n\t"
        "vldd %1, 32(%5)\n\t"
        "vldd %2, 64(%5)\n\t"
        "vldd %3, 96(%5)\n\t"
        "vstd %0,  0(%4)\n\t"
        "vstd %1, 32(%4)\n\t"
        "vstd %2, 64(%4)\n\t"
        "vstd %3, 96(%4)\n\t"
        : "=f"(tmp0), "=f"(tmp1), "=f"(tmp2), "=f"(tmp3)
        : "r"(dest + i), "r"(src + i));
  }
}

void *swmemcpy_vec(void *dest, void *src, size_t cnt){
  for (int i = 0; i < cnt; i += 32){
    intv8 tmp0, tmp1, tmp2, tmp3;
    asm("vldd %0,  0(%2)\n\t"
        "vstd %0,  0(%1)\n\t"
        : "=f"(tmp0)
        : "r"(dest + i), "r"(src + i));
  }
  return dest;
}
extern void *__real_memcpy(void *, void *, size_t);
void *swmemcpy_32u(void *dest, void *src, size_t cnt){
  void *dest_new = dest, *src_new = src;
  int cnt_new = cnt;
  if (cnt_new & 32) {
    intv8 tmp0, tmp1;
    asm("vldw_ul %0,  0(%3)\n\t"
        "vldw_uh %1, 32(%3)\n\t"
        "vbisw %0, %1, %0\n\t"
        "vstw_ul %0,  0(%2)\n\t"
        "vstw_uh %0, 32(%2)\n\t"
        : "=&f"(tmp0), "=&f"(tmp1)
        : "r"(dest_new), "r"(src_new));
    cnt_new -= 32;
    dest_new += 32;
    src_new += 32;
  }
  //return __real_memcpy(dest_new, src_new, cnt_new);
  if (cnt_new & 64){
    intv8 tmp0, tmp1, tmp2, tmp3;
    asm("vldw_ul %0,  0(%5)\n\t"
        "vldw_uh %1, 32(%5)\n\t"
        "vldw_ul %2, 32(%5)\n\t"
        "vldw_uh %3, 64(%5)\n\t"
        "vbisw   %1, %0, %0\n\t"
        "vbisw   %3, %2, %2\n\t"
        "vstw_ul %0,  0(%4)\n\t"
        "vstw_uh %0, 32(%4)\n\t"
        "vstw_ul %2, 32(%4)\n\t"
        "vstw_uh %2, 64(%4)\n\t"
        : "=&f"(tmp0), "=&f"(tmp1), "=&f"(tmp2), "=&f"(tmp3)
        : "r"(dest_new), "r"(src_new));
    cnt_new -= 64;
    dest_new += 64;
    src_new +=64;
  }

  for (int i = 0; i < cnt_new; i += 128){
    intv8 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
    asm("vldw_ul %0,   0(%9)\n\t"
        "vldw_uh %1,  32(%9)\n\t"
        "vldw_ul %2,  32(%9)\n\t"
        "vldw_uh %3,  64(%9)\n\t"
        "vldw_ul %4,  64(%9)\n\t"
        "vldw_uh %5,  96(%9)\n\t"
        "vldw_ul %6,  96(%9)\n\t"
        "vldw_uh %7, 128(%9)\n\t"
        "vbisw   %1, %0, %0\n\t"
        "vbisw   %3, %2, %2\n\t"
        "vbisw   %5, %4, %4\n\t"
        "vbisw   %7, %6, %6\n\t"
        "vstw_ul %0,   0(%8)\n\t"
        "vstw_uh %0,  32(%8)\n\t"
        "vstw_ul %2,  32(%8)\n\t"
        "vstw_uh %2,  64(%8)\n\t"
        "vstw_ul %4,  64(%8)\n\t"
        "vstw_uh %4,  96(%8)\n\t"
        "vstw_ul %6,  96(%8)\n\t"
        "vstw_uh %6, 128(%8)\n\t"
        : "=&f"(tmp0), "=&f"(tmp1), "=&f"(tmp2), "=&f"(tmp3),
          "=&f"(tmp4), "=&f"(tmp5), "=&f"(tmp6), "=&f"(tmp7)
        : "r"(dest_new + i), "r"(src_new + i));
    //simd_loadu(tmp0, src + 1);
    //simd_storeu(tmp0, dest + i);
  }
}

int do_log = 0;
void enable_memcpy_log(){
  do_log = 1;
}
void disable_memcpy_log(){
  do_log = 0;
}
#include <athread.h>
extern void slave_memcpy_cpe(memcpy_arg_t *);
void *__wrap_memcpy(void *dest, void *src, size_t cnt){
  /* if (do_log) */
  /*   printf("%p %p %p %d\n", __builtin_return_address(0), dest, src, cnt); */
  if (cnt < 0x20) return __real_memcpy(dest, src, cnt);
  int align = check_align((long)dest | (long)src | cnt);
  if (cnt > 0x4000 && align >= 4 && athread_idle()) {
    memcpy_arg_t arg;
    arg.dest = dest;
    arg.src = src;
    arg.cnt = cnt;
    athread_spawn(memcpy_cpe, &arg);
    athread_join();
    return dest;
  }

  if (align >= 128) return swmemcpy_vec4(dest, src, cnt);
  if (align >= 32) return swmemcpy_vec(dest, src, cnt);
  if (align >= 4){
    size_t una = cnt & 31L;
    __real_memcpy(dest, src, una);
    swmemcpy_32u(dest + una, src + una, cnt - una);
    return dest;
  }
  return __real_memcpy(dest, src, cnt);
}

#endif
