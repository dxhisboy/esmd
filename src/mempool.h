#ifndef MEMPOOL_H_
#define MEMPOOL_H_
typedef struct mempool {
  void *buffer;
  int *free_list, free_count;
  int block_size;
} mempool_t;
#endif
