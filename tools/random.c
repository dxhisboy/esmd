#include <stdio.h>
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double random(int* idum)
{
  int k;
  double ans;

  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;

  if(*idum < 0) *idum += IM;

  ans = AM * (*idum);
  return ans;
}
int main(){
  int x = 1;
  int i;
  for (i = 0; i <= 5; i ++){
    double r = random(&x);
    printf("%f %d\n", r, x);
  }
}
