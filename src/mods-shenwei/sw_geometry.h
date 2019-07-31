#include <simd.h>
inline doublev4 dsq_atom_box_v4(doublev4 *x_v4, doublev4 *box_o_v4, doublev4 *box_v_v4){
  doublev4 v_v4[3], u_v4[3];
  asm("vsubd %6, %3, %0\n\t"
      "vsubd %7, %4, %1\n\t"
      "vsubd %8, %5, %2\n\t"
      "vcpys $31, %0, %0\n\t"
      "vcpys $31, %1, %1\n\t"
      "vcpys $31, %2, %2\n\t"
      : "=&r"(v_v4[0]), "=&r"(v_v4[0]), "=&r"(v_v4[0])
      : "r"(box_o_v4[0]), "r"(box_o_v4[1]), "r"(box_o_v4[2]),
        "r"(x_v4[0]), "r"(x_v4[1]), "r"(x_v4[2]));
  asm("vsubd %6, %3, %0\n\t"
      "vsubd %7, %4, %1\n\t"
      "vsubd %8, %5, %2\n\t"
      "vsellt %0, $31, %0, %0\n\t"
      "vsellt %1, $31, %1, %1\n\t"
      "vsellt %2, $31, %2, %2\n\t"
      : "=&r"(u_v4[0]), "=&r"(u_v4[0]), "=&r"(u_v4[0])
      : "r"(box_v_v4[0]), "r"(box_v_v4[1]), "r"(box_v_v4[2]),
        "r"(v_v4[0]), "r"(v_v4[1]), "r"(v_v4[2]));
  return u_v4[0] * u_v4[0] + u_v4[1] * u_v4[1] + u_v4[2] * u_v4[2];
}
