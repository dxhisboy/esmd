#include <memory.h>
#include <multiproc.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
//#define DEBUG_THIS_FILE
#include <log.h>
inline long current_usec(){
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * 1000000 + tv.tv_usec;
}
typedef struct field {
  long sum, min, max;
  int argmin, argmax, argn;
} field_t;
typedef struct timerec {
  char name[256];
  field_t time, cnt;
} timerec_t;
#define PPH_ONETIME
#define PPH_NAME timertab
#define PPH_TYPE timerec_t
#define PPH_HASH(x) (pph_hash_string((x)->name))
#define PPH_EQ(x, y) (!strcmp((x)->name, (y)->name))
#define PPH_CALLOC(x, y) esmd_calloc(x, y, "counter tab")
#define PPH_ON_DELETE(x) esmd_free(x)
#define PPH_FREE esmd_free

#include <pphash.h>
static timertab_t timertab;
void timer_init(){
  timertab_init(&timertab, 0);
}

void timer_start(const char *tname) {
  timerec_t *cast_name = (void*)tname;
  timerec_t **slot = timertab_find_slot(&timertab, cast_name);
  assert(slot);
  timerec_t *timer;
  if (pph_slot_is_empty(slot)){
    timer = esmd_malloc(sizeof(timerec_t), "timertab entry");
    strcpy(timer->name, tname);
    timer->time.sum = 0;
    timer->cnt.sum = 0;
    timertab_insert(&timertab, slot, timer);
  } else {
    timer = *slot;
  }
  timer->time.sum -= current_usec();
  timer->cnt.sum ++;
}

void timer_stop(const char *tname) {
  timerec_t *cast_name = (void*)tname;
  timerec_t **slot = timertab_find_slot(&timertab, cast_name);
  assert(slot);
  timerec_t *timer;
  if (pph_slot_is_empty(slot)){
    timer = esmd_malloc(sizeof(timerec_t), "timertab entry");
    strcpy(timer->name, tname);
    timer->time.sum = 0;
    timertab_insert(&timertab, slot, timer);
  } else {
    timer = *slot;
  }
  timer->time.sum += current_usec();
}

static void maintain_timer(timerec_t *timer, void *arg) {
  int rank = (int)(long)arg;
  timer->time.min = timer->time.sum;
  timer->time.max = timer->time.sum;
  timer->time.argmin = rank;
  timer->time.argmax = rank;
  timer->time.argn = 1;
  
  timer->cnt.min = timer->cnt.sum;
  timer->cnt.max = timer->cnt.sum;
  timer->cnt.argmin = rank;
  timer->cnt.argmax = rank;
  timer->cnt.argn = 1;
}

static void merge_timer(timerec_t *timer, void *gbl) {
  timertab_t *gbl_tab = (timertab_t*)gbl;
  timerec_t **slot = timertab_find_slot(gbl_tab, timer);
  assert(slot);
  if (pph_slot_is_empty(slot)){
    timerec_t *gbl_timer = esmd_malloc(sizeof(timerec_t), "gbltimertab entry");
    memcpy(gbl_timer, timer, sizeof(timerec_t));
    timertab_insert(gbl_tab, slot, gbl_timer);
  } else {
    timerec_t *gbl_timer = *slot;
    gbl_timer->time.sum += timer->time.sum;
    gbl_timer->time.argn += timer->time.argn;
    if (timer->time.min < gbl_timer->time.min){
      gbl_timer->time.min = timer->time.min;
      gbl_timer->time.argmin = timer->time.argmin;
    }
    if (timer->time.max > gbl_timer->time.max){
      gbl_timer->time.max = timer->time.max;
      gbl_timer->time.argmax = timer->time.argmax;
    }
    gbl_timer->cnt.sum += timer->cnt.sum;
    gbl_timer->cnt.argn += timer->cnt.argn;
    if (timer->cnt.min < gbl_timer->cnt.min){
      gbl_timer->cnt.min = timer->cnt.min;
      gbl_timer->cnt.argmin = timer->cnt.argmin;
    }
    if (timer->cnt.max > gbl_timer->cnt.max){
      gbl_timer->cnt.max = timer->cnt.max;
      gbl_timer->cnt.argmax = timer->cnt.argmax;
    }

  }
}

#define allreduce_tag 400
static void timer_allreduce(timertab_t *gbl_timertab, MPI_Comm comm){
  int np, me;
  MPI_Comm_size(comm, &np);
  MPI_Comm_rank(comm, &me);
  timertab_init(gbl_timertab, 0);
  timertab_traverse(&timertab, maintain_timer, (void*)(long)me);
  timertab_traverse(&timertab, merge_timer, gbl_timertab);

  int n = gbl_timertab->cnt, gbl_n;
  MPI_Allreduce(&n, &gbl_n, 1, MPI_INT, MPI_MAX, comm);
  timerec_t *ent_recv = esmd_malloc(sizeof(timerec_t) * gbl_n, "perf reduce ent");
  timerec_t *ent_send = esmd_malloc(sizeof(timerec_t) * gbl_n, "perf reduce ent");
  for (int stride = 1; stride < np; stride += stride) {
    if (me & stride){
      debug("%d send-\n", me);
      int nrec = timertab_pack(gbl_timertab, ent_send);
      MPI_Send(ent_send, sizeof(timerec_t) * nrec, MPI_CHAR, me - stride, 400, comm);
      debug("%d send\n", me);
      break;
    } else if (me + stride < np){
      debug("%d recv-\n", me);
      MPI_Status stat;
      MPI_Recv(ent_recv, sizeof(timerec_t) * gbl_n, MPI_CHAR, me + stride, 400, comm, &stat);
      int recv_count;
      MPI_Get_count(&stat, MPI_CHAR, &recv_count);
      int nrec = recv_count / sizeof(timerec_t);
      for (int i = 0; i < nrec; i ++){
	merge_timer(ent_recv + i, gbl_timertab);
      }
      /* if (me == 0){ */
      /* 	timerec_t **slot = timertab_find_slot(gbl_timertab, (const timerec_t*)("force")); */
      /* 	debug("%d\n", (*slot)->cnt.sum); */
      /* } */
      debug("%d recv\n", me);
    }
  }
  esmd_free(ent_recv);
  esmd_free(ent_send);
}

inline void timerec_print(timerec_t *timer, void *arg){
  double avg = timer->time.sum * 1e-6 / timer->time.argn;
  long cnt = timer->cnt.sum;
  double max = timer->time.max * 1e-6;
  double min = timer->time.min * 1e-6;
  int argmax = timer->time.argmax;
  int argmin = timer->time.argmin;
  printf("%-20s%10ld%10.3f%10.3f(%6d)%10.3f(%6d)\n",
	 timer->name, cnt, avg,
	 max, argmax, min, argmin);
}
void timer_print(MPI_Comm comm){
  timertab_t gbl_timertab;
  timer_allreduce(&gbl_timertab, comm);
  int me;
  MPI_Comm_rank(comm, &me);
  if (me == 0){
    printf("%-20s%10s%10s%10s(%6s)%10s(%6s)\n",
	   "timer", "count", "mean",
	   "max", "rank", "min", "rank");

    timertab_traverse(&gbl_timertab, timerec_print, NULL);
  }
  timertab_destroy(&gbl_timertab);
}
