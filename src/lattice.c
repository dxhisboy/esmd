#include <math.h>

#include <lattice.h>
#define IA   16807
#define IM   2147483647
#define AM   (1.0/IM)
#define IQ   127773
#define IR   2836
#define MASK 123459876

double random(int *seed){
  int k = *seed / IQ;
  *seed = IA * (*seed - k * IQ) - IR * k;
  if (*seed < 0){
    *seed += IM;
  }
  return AM * *seed;
}

areal offset_fcc[4][3] = {
  {0.0, 0.0, 0.0},
  {0.5, 0.5, 0.0},
  {0.5, 0.0, 0.5},
  {0.0, 0.5, 0.5}
};

int weight_fcc = 4;

lattice_t lattices[] = {
  //fcc
  {
    {
      {0.0, 0.0, 0.0},
      {0.5, 0.5, 0.0},
      {0.5, 0.0, 0.5},
      {0.0, 0.5, 0.5}
    }, 4
  }
};
static double lattice_length(esmd_t *md, lattice_type type, unit_type utype){
  if (utype == UNIT_LJ){
    int weight = lattices[type].weight;
    double scale = pow(weight / dens, 1 / 3.0);

  }
}
