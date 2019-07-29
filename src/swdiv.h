#ifndef SWDIV_H_
#define SWDIV_H_
typedef struct {
  int64_t k, m;
} div_idx_t;
inline void build_div_idx(div_idx_t *idx, uint32_t c){
  for (int k = 0; k < 32; k ++){
    if ((1L << k) <= c && (1L << k + 1) > c) {
      idx->k = k;
      break;
    }
  }
  if ((1 << idx->k) == c) {
    idx->m = 1;
  } else {
    idx->k = idx->k + 32;
    idx->m = ((1L << idx->k) + c - 1) / c;
  }
}
inline uint32_t div_idx(uint32_t x, div_idx_t *idx){
  return (x * idx->m) >> idx->k;
}

#endif
