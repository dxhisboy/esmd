#include <data.h>
enum lattice_type {
  LAT_FCC
};
enum unit_type {
  UNIT_LJ,
  UNIT_REAL
};

typedef struct lattice {
  areal (*offset)[3];
  areal lx, ly, lz;
  int weight;
} lattice_t;
