#ifndef TIMER_H_
#define TIMER_H_
#include <data.h>
//function signatures
void timer_init();
void timer_start(const char *tname);
void timer_stop(const char *tname);
void timer_print(MPI_Comm comm);
//end function signatures
#endif
