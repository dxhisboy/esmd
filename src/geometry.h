#include <math.h>
#include <util.h>
inline areal dsq_atom_box(areal *x, areal *box_o, areal *box_v){
  areal v[3], u[3];
  v[0] = fabs(x[0] - box_o[0]);
  v[1] = fabs(x[1] - box_o[1]);
  v[2] = fabs(x[2] - box_o[2]);
  u[0] = max(v[0] - box_v[0], 0);
  u[1] = max(v[1] - box_v[1], 0);
  u[2] = max(v[2] - box_v[2], 0);
  return u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
}
