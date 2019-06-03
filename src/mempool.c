#include <stdlib.h>
#include <mempool.h>
void mempool_init(mempool_t *pool, int block_size, int num_blocks) {
  pool->block_size = block_size;
  pool->buffer = malloc(block_size * num_blocks);
  pool->free_list = (int*)malloc(num_blocks * sizeof(int));
  pool->free_count = num_blocks;

  for (int i = 0; i < num_blocks; i ++){
    pool->free_list[i] = i * block_size;
  }
}
void *mempool_get(mempool_t *pool) {
  int block_head = pool->free_list[-- pool->free_count];
  return pool->buffer + block_head;
}

void mempool_return(mempool_t *pool, void *ptr) {
  int block_head = ptr - pool->buffer;
  pool->free_list[pool->free_count ++] = block_head;
}

void mempool_destroy(mempool_t *pool) {
  free(pool->free_list);
  free(pool->buffer);
}
