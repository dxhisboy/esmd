#ifndef MEMORY_H_
#define MEMORY_H_
#include <stdlib.h>
#include <cppdefs.h>
typedef struct mempool {
  void *buffer;
  int *free_list, free_count;
  int block_size;
} mempool_t;

void mempool_init(mempool_t*, int, int, char *);
void *mempool_get(mempool_t*);
void mempool_return(mempool_t*, void *);
void mempool_destroy(mempool_t*);
void *esmd_malloc(size_t, char*);
void esmd_free(void*);
void memory_init();
void memory_print();
#endif
