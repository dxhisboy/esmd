#ifndef PPHASH_SHARED
#define PPHASH_SHARED
#define HTAB_EMPTY ((void*)0)
#define HTAB_DELETED ((void*)1)
typedef unsigned int hashval_t;
struct prime_ent
{
  hashval_t prime;
  hashval_t inv;
  hashval_t inv_m2;	/* inverse of prime-2 */
  hashval_t shift;
};

static struct prime_ent const prime_tab[] = {
  {          7, 0x24924925, 0x9999999b, 2 },
  {         13, 0x3b13b13c, 0x745d1747, 3 },
  {         31, 0x08421085, 0x1a7b9612, 4 },
  {         61, 0x0c9714fc, 0x15b1e5f8, 5 },
  {        127, 0x02040811, 0x0624dd30, 6 },
  {        251, 0x05197f7e, 0x073260a5, 7 },
  {        509, 0x01824366, 0x02864fc8, 8 },
  {       1021, 0x00c0906d, 0x014191f7, 9 },
  {       2039, 0x0121456f, 0x0161e69e, 10 },
  {       4093, 0x00300902, 0x00501908, 11 },
  {       8191, 0x00080041, 0x00180241, 12 },
  {      16381, 0x000c0091, 0x00140191, 13 },
  {      32749, 0x002605a5, 0x002a06e6, 14 },
  {      65521, 0x000f00e2, 0x00110122, 15 },
  {     131071, 0x00008001, 0x00018003, 16 },
  {     262139, 0x00014002, 0x0001c004, 17 },
  {     524287, 0x00002001, 0x00006001, 18 },
  {    1048573, 0x00003001, 0x00005001, 19 },
  {    2097143, 0x00004801, 0x00005801, 20 },
  {    4194301, 0x00000c01, 0x00001401, 21 },
  {    8388593, 0x00001e01, 0x00002201, 22 },
  {   16777213, 0x00000301, 0x00000501, 23 },
  {   33554393, 0x00001381, 0x00001481, 24 },
  {   67108859, 0x00000141, 0x000001c1, 25 },
  {  134217689, 0x000004e1, 0x00000521, 26 },
  {  268435399, 0x00000391, 0x000003b1, 27 },
  {  536870909, 0x00000019, 0x00000029, 28 },
  { 1073741789, 0x0000008d, 0x00000095, 29 },
  { 2147483647, 0x00000003, 0x00000007, 30 },
  /* Avoid "decimal constant so large it is unsigned" for 4294967291.  */
  { 0xfffffffb, 0x00000006, 0x00000008, 31 }
};

static inline hashval_t
pph_hash_string(const char *s)
{
  hashval_t ret = 0;
  const char *p;
  for (p = s; *p; p ++){
    ret = ret * 67 + (*p) - 113;
  }
  return ret;
}

#define mix(a,b,c)                                      \
  {                                                     \
    a -= b; a -= c; a ^= (c>>13);                       \
    b -= c; b -= a; b ^= (a<< 8);                       \
    c -= a; c -= b; c ^= ((b&0xffffffff)>>13);          \
    a -= b; a -= c; a ^= ((c&0xffffffff)>>12);          \
    b -= c; b -= a; b = (b ^ (a<<16)) & 0xffffffff;     \
    c -= a; c -= b; c = (c ^ (b>> 5)) & 0xffffffff;     \
    a -= b; a -= c; a = (a ^ (c>> 3)) & 0xffffffff;     \
    b -= c; b -= a; b = (b ^ (a<<10)) & 0xffffffff;     \
    c -= a; c -= b; c = (c ^ (b>>15)) & 0xffffffff;     \
  }

typedef __intptr_t intptr_t;
static hashval_t
pph_hash_pointer (const void *p)
{
  intptr_t v = (intptr_t) p;
  unsigned a, b, c;
  a = b = 0x9e3779b9;
  a += v >> (sizeof (intptr_t) * 8 / 2);
  b += v & (((intptr_t) 1 << (sizeof (intptr_t) * 8 / 2)) - 1);
  c = 0x42135234;
  mix (a, b, c);
  return c;
}

static inline int
pph_slot_is_empty(const void* slot){
  const void *entry = *((void**)slot);
  return entry == HTAB_EMPTY || entry == HTAB_DELETED;
}

#endif

#ifndef PPH_CALLOC
#define PPH_CALLOC calloc
#endif
#ifndef PPH_FREE
#define PPH_FREE free
#endif

#ifndef PPH_ON_DELETE
#define PPH_ON_DELETE(x)
#endif

#ifndef PPH_ENTCPY
#define PPH_ENTCPY(x, y) memcpy(x, y, sizeof(PPH_TYPE));
#endif

#define __CAT__(x, y) x ## _ ## y
#define CAT(x, y) __CAT__(x, y)

struct PPH_NAME {
  PPH_TYPE **slots;
  int cap, cnt;
  long nquery, nscan;
  const struct prime_ent *prime_cur;
};
typedef struct PPH_NAME CAT(PPH_NAME, t);
#define htab_t struct PPH_NAME
//struct prime_ent p = {127, 0x02040811, 0x0624dd30, 6};
static inline void CAT(PPH_NAME, check_cap)(htab_t *htab);
static inline PPH_TYPE**
CAT(PPH_NAME, find_slot_hash)(htab_t *htab, const PPH_TYPE *element, hashval_t hash)
{
  CAT(PPH_NAME, check_cap)(htab);
  htab->nquery ++;

  hashval_t cap = htab->cap;
  PPH_TYPE **avail = NULL;

  hashval_t index = hash % htab->prime_cur->prime;
  PPH_TYPE **slot = htab->slots + index;

  if (*slot == HTAB_EMPTY)
    return slot;
  else if (*slot == HTAB_DELETED)
    avail = slot;
  else if (PPH_EQ(*slot, element))
    return slot;
  
  hashval_t hash2 = 1 + hash % (htab->prime_cur->prime - 2);
  while (1){
    htab->nscan ++;
    index += hash2;
    if (index >= cap)
      index -= cap;
    slot = htab->slots + index;
    if (*slot == HTAB_EMPTY) 
      return avail ? avail : slot;
    else if (*slot == HTAB_DELETED)
      avail = avail ? avail : slot;
    else if (PPH_EQ(*slot, element))
      return slot;
  }
}

static inline PPH_TYPE**
CAT(PPH_NAME, find_slot)(htab_t *htab, const PPH_TYPE *element){
  hashval_t hv = PPH_HASH(element);
  return CAT(PPH_NAME, find_slot_hash)(htab, element, hv);
}
static inline PPH_TYPE*
CAT(PPH_NAME, find_entry)(htab_t *htab, const PPH_TYPE *element, hashval_t hash)
{
  htab->nquery ++;
  hashval_t cap = htab->cap;
  PPH_TYPE **avail = NULL;

  hashval_t index = hash % htab->prime_cur->prime;
  PPH_TYPE *entry = htab->slots[index];

  entry = htab->slots[index];
  if (entry == HTAB_EMPTY || entry != HTAB_DELETED && PPH_EQ(entry, element))
    return entry;
  
  hashval_t hash2 = 1 + hash % (htab->prime_cur->prime - 2);
  
  while (1){
    htab->nscan ++;
    index += hash2;
    if (index >= cap)
      index -= cap;
    entry = htab->slots[index];
    if (entry == HTAB_EMPTY || entry != HTAB_DELETED && PPH_EQ(entry, element))
      return entry;
  }
}

static inline PPH_TYPE*
CAT(PPH_NAME, init)(htab_t *htab, int initial_size)
{
  htab->prime_cur = prime_tab;
  while (htab->prime_cur->prime * 3 < initial_size * 4)
    htab->prime_cur ++;
  
  htab->nquery = 0;
  htab->nscan = 0;
  htab->cap = htab->prime_cur->prime;
  htab->cnt = 0;
  htab->slots = PPH_CALLOC(htab->cap, sizeof(PPH_TYPE*));
}

static inline void
CAT(PPH_NAME, traverse)(htab_t *htab, void (*callback)(PPH_TYPE*, void *), void *arg)
{
  PPH_TYPE **top = htab->slots + htab->cap, **slot;
  for (slot = htab->slots; slot != top; slot ++){
    if (*slot != HTAB_EMPTY && *slot != HTAB_DELETED) {
      callback(*slot, arg);
    }
  }
}

static inline void
CAT(PPH_NAME, check_cap)(htab_t *htab){
  if (htab->cap * 3 < htab->cnt * 4){
    puts("extending");
    htab->prime_cur ++;
    PPH_TYPE **old_slots = htab->slots;
    PPH_TYPE **top = htab->slots + htab->cap, **slot;
    htab->cap = htab->prime_cur->prime;
    htab->slots = PPH_CALLOC(htab->cap, sizeof(PPH_TYPE*));
    int nquery_old = htab->nquery;
    int nscan_old = htab->nscan;
    for (slot = old_slots; slot != top; slot ++){
      if (*slot != HTAB_EMPTY && *slot != HTAB_DELETED){
	hashval_t hv = PPH_HASH(*slot);
	PPH_TYPE **new_slot = CAT(PPH_NAME, find_slot_hash)(htab, *slot, hv);
	*new_slot = *slot;
      }
    }
    PPH_FREE(old_slots);
    htab->nquery = nquery_old;
    htab->nscan = nscan_old;
  }
}

static inline void
CAT(PPH_NAME, insert)(htab_t *htab, PPH_TYPE **slot, PPH_TYPE *element){
  if (*slot == HTAB_EMPTY) htab->cnt ++;
  *slot = element;
}

static inline void
CAT(PPH_NAME, delete)(htab_t *htab, PPH_TYPE **slot){
  PPH_ON_DELETE(*slot);
  *slot = HTAB_DELETED;
}

static inline int
CAT(PPH_NAME, pack)(htab_t *htab, PPH_TYPE *entries){
  int nrec_write = 0;
  PPH_TYPE **top = htab->slots + htab->cap, **slot;
  for (slot = htab->slots; slot != top; slot ++){
    if (*slot != HTAB_EMPTY && *slot != HTAB_DELETED) {
      PPH_ENTCPY(entries + nrec_write, *slot);
      nrec_write ++;
    }
  }
  return nrec_write;
}

static inline void
CAT(PPH_NAME, destroy)(htab_t *htab){
  PPH_TYPE **top = htab->slots + htab->cap, **slot;
  for (slot = htab->slots; slot != top; slot ++){
    if (*slot != HTAB_EMPTY && *slot != HTAB_DELETED) {
      PPH_ON_DELETE(*slot);
    }
  }
  PPH_FREE(htab->slots);
}
#undef CAT
#undef __CAT__
#undef htab_t
#ifdef PPH_ONETIME
#undef PPH_NAME
#undef PPH_HASH
#undef PPH_TYPE
#undef PPH_EQ
#undef PPH_CALLOC
#undef PPH_FREE
#undef PPH_ON_DELETE
#undef PPH_ENTCPY
#endif
