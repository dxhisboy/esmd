#ifndef MEMORY_H_
#define MEMORY_H_
#include <stdlib.h>
typedef struct mempool {
  void *buffer;
  int *free_list, free_count;
  int block_size;
} mempool_t;

//function signatures
void mempool_init(mempool_t *pool, int block_size, int num_blocks, const char *name);
void *mempool_get(mempool_t *pool);
void mempool_return(mempool_t *pool, void *ptr);
void mempool_destroy(mempool_t *pool);
void memory_init();
void memory_print();
void *esmd_malloc(size_t size, const char *name);
void *esmd_calloc(size_t nmemb, size_t size, char *name);
void esmd_free(void *ptr);
//end function signatures
#endif
