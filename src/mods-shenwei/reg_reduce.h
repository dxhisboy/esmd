#ifndef REG_REDUCE_H_
#define REG_REDUCE_H_
#include <simd.h>
inline doublev4 reg_reduce_doublev4(doublev4 in){
  doublev4 myresult = in;
  for (int stride = 1; stride < 8; stride += stride){
    if (_COL & stride){
      simd_putr(myresult, _COL - stride);
      break;
    } else {
      doublev4 tmp;
      tmp = simd_getr(tmp);
      myresult += tmp;
    }
    //athread_syn(ROW_SCOPE, 0xff);
  }
  if (_COL == 0){
    for (int stride = 1; stride < 8; stride += stride){
      if (_ROW & stride){
        simd_putc(myresult, _ROW - stride);
        break;
      } else {
        doublev4 tmp;
        tmp = simd_getc(tmp);
        myresult += tmp;
      }
    }
  }
  return myresult;
}
#endif
