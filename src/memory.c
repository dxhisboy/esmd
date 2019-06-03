#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <hashtab.h>

#include <data.h>
typedef struct memrec {
  char name[256];
  size_t size, cnt;
} memrec_t;

typedef struct meminfo {
  memrec_t *rec;
  void *raw;
  size_t size;
} meminfo_t;

static htab_t htab;
memrec_t htab_mr = {"memory hashtab", 0, 0};
mempool_t rec_pool;
extern void *camd_malloc(size_t, char*);
void mempool_init(mempool_t *pool, int block_size, int num_blocks, char *name) {
  pool->block_size = block_size;
  pool->buffer = camd_malloc(block_size * num_blocks, name);
  pool->free_list = camd_malloc(sizeof(int) * num_blocks, name);
  pool->free_count = num_blocks;

  for (int i = 0; i < num_blocks; i ++){
    pool->free_list[i] = i * block_size;
  }
}

static void mempool_init_prealloc(mempool_t *pool, int block_size, int num_blocks, void *buffer, int *list){
  pool->block_size = block_size;
  pool->buffer = buffer;
  pool->free_list = list;
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
  camd_free(pool->free_list);
  camd_free(pool->buffer);
}

static int mr_eq(const void *e1, const void *e2) {
  const char *r1 = e1, *r2 = e2;
  return strcmp(r1, r2) == 0;
}

static hashval_t mr_hash(const void *e) {
  const char *r = e;
  return htab_hash_string(r);
}

static void mr_del(void *e){
}

static int mr_print_trav(void **mrpp, void *help){
  memrec_t *node = *mrpp;
  printf("%s: %d bytes malloced in %d mallocs\n", node->name, node->size, node->cnt);
  return 1;
}

static void *camd_aligned_malloc(size_t size){
  long raw = (long)malloc(size + sizeof(meminfo_t) + MEMORY_ALIGN_MASK);
  void *ret = (void*)((raw + sizeof(meminfo_t) + MEMORY_ALIGN_MASK) & ~MEMORY_ALIGN_MASK);
  meminfo_t *info = ret - sizeof(meminfo_t);
  info->raw = (void*)raw;
  info->size = size;
  info->rec = NULL;
  return ret;
}

static void camd_aligned_free(void *ptr){
  meminfo_t *info = ptr - sizeof(meminfo_t);
  free(info->raw);
}
static void *camd_htab_calloc(size_t size, size_t count){
  htab_mr.size += size * count;
  htab_mr.cnt += 1;
  void *ret = camd_aligned_malloc(size * count);
  memset(ret, 0, size * count);
  return ret;
}

static void camd_htab_free(void *ptr){
  meminfo_t *info = ptr - sizeof(meminfo_t);
  htab_mr.size -= info->size;
  htab_mr.cnt --;
  camd_aligned_free(ptr);
}

void memory_init(){
  htab = htab_create_alloc(N_MEMREC, mr_hash, mr_eq, mr_del, camd_htab_calloc, camd_htab_free);
  void **slot = htab_find_slot(htab, &htab_mr, INSERT);
  assert(*slot == NULL);
  *slot = &htab_mr;

  void *rec_buffer = camd_aligned_malloc(sizeof(memrec_t) * N_MEMREC);
  meminfo_t *buffer_info = rec_buffer - sizeof(meminfo_t);
  buffer_info->rec = &htab_mr;
  
  void *rec_list = camd_aligned_malloc(sizeof(int) * N_MEMREC);
  meminfo_t *list_info = rec_list - sizeof(meminfo_t);
  list_info->rec = &htab_mr;
  
  mempool_init_prealloc(&(rec_pool), sizeof(memrec_t), N_MEMREC, rec_buffer, rec_list);
  htab_mr.size += sizeof(int) * N_MEMREC + sizeof(memrec_t) * N_MEMREC;
  htab_mr.cnt += 1;
}

void memory_print(){
  htab_traverse(htab, mr_print_trav, NULL);
}

void *camd_malloc(size_t size, char *name){
  void *ret = camd_aligned_malloc(size);
  meminfo_t *info = ret - sizeof(meminfo_t);
  void **slot = htab_find_slot(htab, name, INSERT);
  assert(slot);
  if (*slot == NULL){
    memrec_t *rec = mempool_get(&(rec_pool));
    strcpy(rec->name, name);
    rec->size = size;
    rec->cnt = 1;
    *slot = rec;
    info->rec = rec;
  } else {
    memrec_t *rec = *slot;
    rec->size += size;
    rec->cnt += 1;
    info->rec = rec;
  }
}

void camd_free(void *ptr){
  meminfo_t *info = ptr - sizeof(meminfo_t);
  info->rec->size -= info->size;
  info->rec->cnt --;
  camd_aligned_free(ptr);
}

#ifdef TEST
int main(){
  memory_init();
  //mempool_destroy(&(rec_pool));
  memory_print();
}
#endif
