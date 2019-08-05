#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
//#include <hashtab.h>

#include <data.h>
//#define DEBUG_THIS_FILE
#include <log.h>
typedef struct memrec {
  char name[256];
  size_t size, cnt;
} memrec_t;

typedef struct meminfo {
  memrec_t *rec;
  void *raw;
  size_t size;
} meminfo_t;

memrec_t mr_self = {"memory hashtab", 0, 0};
mempool_t rec_pool;

static void *esmd_htab_calloc(size_t count, size_t size);
static void esmd_htab_free(void *ptr);

#define PPH_NAME mrtab
#define PPH_TYPE memrec_t
#define PPH_HASH(x) (pph_hash_string((x)->name))
#define PPH_EQ(x, y) (!strcmp((x)->name, (y)->name))
#define PPH_CALLOC esmd_htab_calloc
#define PPH_FREE esmd_htab_free
#include "pphash.h"

//static htab_t htab;
static mrtab_t mrtab;
extern void *esmd_malloc(size_t, const char*);
void mempool_init(mempool_t *pool, int block_size, int num_blocks, const char *name) {
  pool->block_size = block_size;
  pool->buffer = esmd_malloc(block_size * num_blocks, name);
  pool->free_list = esmd_malloc(sizeof(int) * num_blocks, name);
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
  esmd_free(pool->free_list);
  esmd_free(pool->buffer);
}

static void mr_print_trav(memrec_t *node, void *help){
  char *units[] = {"B", "K", "M", "G", "T", "P", "E"};
  int unit_ptr = 0;
  size_t size = node->size;
  while (size > 1024){
    size /= 1024;
    unit_ptr ++;
  }
  printf("%4ld%s %6ld %s\n", size, units[unit_ptr], node->cnt, node->name);
}

static void *esmd_aligned_malloc(size_t size){
  long raw = (long)malloc(size + sizeof(meminfo_t) + MEMORY_ALIGN_MASK);
  void *ret = (void*)((raw + sizeof(meminfo_t) + MEMORY_ALIGN_MASK) & ~MEMORY_ALIGN_MASK);
  if (raw == 0){
    perror("Malloc failed");
    exit(1);
  }
  debug("%p\n", ret);
  meminfo_t *info = ret - sizeof(meminfo_t);
  info->raw = (void*)raw;
  info->size = size;
  info->rec = NULL;
  return ret;
}

static void esmd_aligned_free(void *ptr) {
  meminfo_t *info = ptr - sizeof(meminfo_t);
  free(info->raw);
}

static void *esmd_htab_calloc(size_t count, size_t size) {
  mr_self.size += size * count;
  mr_self.cnt += 1;
  void *ret = esmd_aligned_malloc(size * count);
  memset(ret, 0, size * count);
  return ret;
}

static void esmd_htab_free(void *ptr) {
  meminfo_t *info = ptr - sizeof(meminfo_t);
  mr_self.size -= info->size;
  mr_self.cnt --;
  esmd_aligned_free(ptr);
}

void memory_init() {
  mrtab_init(&mrtab, 0);
  memrec_t **slot = mrtab_find_slot(&mrtab, &mr_self);
  assert(*slot == NULL);
  *slot = &mr_self;

  void *rec_buffer = esmd_aligned_malloc(sizeof(memrec_t) * N_MEMREC);
  meminfo_t *buffer_info = rec_buffer - sizeof(meminfo_t);
  buffer_info->rec = &mr_self;
  
  void *rec_list = esmd_aligned_malloc(sizeof(int) * N_MEMREC);
  meminfo_t *list_info = rec_list - sizeof(meminfo_t);
  list_info->rec = &mr_self;
  
  mempool_init_prealloc(&(rec_pool), sizeof(memrec_t), N_MEMREC, rec_buffer, rec_list);

  mr_self.size += sizeof(int) * N_MEMREC + sizeof(memrec_t) * N_MEMREC;
  mr_self.cnt += 1;
}

void memory_print() {
  mrtab_traverse(&mrtab, mr_print_trav, NULL);
}

void *esmd_malloc(size_t size, const char *name) {
  void *ret = esmd_aligned_malloc(size);
  meminfo_t *info = ret - sizeof(meminfo_t);
  //void **slot = htab_find_slot(htab, name, INSERT);
  memrec_t *cast_name = (void*)name;
  memrec_t **slot = mrtab_find_slot(&mrtab, cast_name);
  assert(slot);
  if (pph_slot_is_empty(slot)){
    memrec_t *rec = mempool_get(&(rec_pool));
    strcpy(rec->name, name);
    rec->size = size;
    rec->cnt = 1;
    //*slot = rec;
    mrtab_insert(&mrtab, slot, rec);
    info->rec = rec;
  } else {
    memrec_t *rec = *slot;
    rec->size += size;
    rec->cnt += 1;
    info->rec = rec;
  }
  return ret;
}

void *esmd_calloc(size_t nmemb, size_t size, char *name) {
  void *ret = esmd_malloc(nmemb * size, name);
  memset(ret, 0, nmemb * size);
  return ret;
}
void esmd_free(void *ptr) {
  meminfo_t *info = ptr - sizeof(meminfo_t);
  info->rec->size -= info->size;
  info->rec->cnt --;
  esmd_aligned_free(ptr);
}

