#include <memory.h>
#include <multiproc.h>
enum counters{
  CNT_USEC,
  N_CNTR
};
#include <sys/time.h>
void read_counters(long *dest){
  struct timeval tv;
  gettimeofday(&tv);
  dest[CNT_USEC] = tv.tv_sec * 1000000 + tv.tv_usec;
}
typedef struct evtent {
  char name[256];
  long evts[N_CNTR];
} evtent_t;
#define PPH_NAME evttab
#define PPH_TYPE counter_t
#define PPH_HASH(x) (pph_hash_string((x)->name))
#define PPH_EQ(x, y) (!strcmp((x)->name, (y)->name))
#define PPH_CALLOC(x, y) esmd_calloc(x, y, "counter tab")
#define PPH_FREE esmd_free

#include <pphash.h>
static evttab_t evttab;
void perf_init(int pid){
  evttab_init(&evttab, 0);
}

void perf_start(const char *tname) {
  evtent_t *cast_name = (void*)tname;
  evtent_t **slot = evttab_find_slot(&evttab, cast_name);
  assert(slot);
  evtent_t *evt;
  if (pph_slot_is_empty(slot)){
    evt = esmd_malloc(sizeof(evtent_t));
    strcpy(evt->name, tname);
    evt->sum = 0;
    for (int i = 0; i < N_EVTR; i ++)
      evt->evts[i] = 0;
  } else {
    evt = *slot;
  }
  long evts[N_EVTR];
  read_counters(evts);
  for (int i = 0; i < N_EVTR; i ++)
    evt.evt[i] -= evts[i];
}

void perf_stop(const char *tname) {
  evtent_t *cast_name = (void*)tname;
  evtent_t **slot = evttab_find_slot(&evttab, cast_name);
  assert(slot);
  evtent_t *evt;
  if (pph_slot_is_empty(slot)){
    evt = esmd_malloc(sizeof(evtent_t));
    strcpy(evt->name, tname);
    evt->sum = 0;
    evttab_insert(&evttab, slot, evt);
  } else {
    evt = *slot;
  }
  long evts[N_EVTR];
  read_counters(evts);
  for (int i = 0; i < N_EVTR; i ++)
    evt.evt[i] += evts[i];
}

void perf_allreduce(){
  
}
