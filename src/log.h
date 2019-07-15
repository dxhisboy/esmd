#ifdef DEBUG_THIS_FILE
#define debug(fmt, ...) fprintf(stderr, "\033[1;36m[DBG]\033[0m " fmt, __VA_ARGS__)
#else
#define debug(fmt, ...)
#endif

#define info(fmt, ...) fprintf(stdout, "\033[1;32m[INF]\033[0m " fmt, __VA_ARGS__)
#define error(fmt, ...) fprintf(stderr, "\033[1;31m[ERR]\033[0m " fmt, __VA_ARGS__)

