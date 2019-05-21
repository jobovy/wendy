/*
  parallel_sort.c: A parallel version of sort
 */
#ifndef __PARALLEL_SORT_H__
#define __PARALLEL_SORT_H__
#include<stdlib.h>
#ifndef PARALLEL_SERIAL_SORT_SWITCH
#define PARALLEL_SERIAL_SORT_SWITCH 10000
#endif
#ifndef PARALLEL_SERIAL_MERGE_SWITCH
#define PARALLEL_SERIAL_MERGE_SWITCH 50000
#endif
// Functions
void parallel_sort(void *,size_t,size_t,
		   int (*)(const void *, const void *));
#endif /* parallel_sort.h */
