/*
  parallel_sort.c: A parallel version of sort, not stable
 */
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include <parallel_sort.h>

void serial_merge(void *src1, void *src2,
		  void *dest,
		  size_t nmemb1, size_t nmemb2, size_t size, 
		  int (*compar ) (const void *, const void * )){
  register void* left_ptr= src1;
  register void* left_ptr_max= src1+size*nmemb1;
  register void* right_ptr= src2;
  register void* right_ptr_max= src2+size*nmemb2;
  // Special case where the arrays don't overlap
  if (compar(left_ptr_max-size,right_ptr) < 0) {
    memcpy(dest,left_ptr,left_ptr_max-left_ptr);
    dest+= left_ptr_max-left_ptr;
    memcpy(dest,right_ptr,right_ptr_max-right_ptr);
    return;
  } else if (compar(right_ptr_max-size,left_ptr) < 0) {
    memcpy(dest,right_ptr,right_ptr_max-right_ptr);
    dest+= right_ptr_max-right_ptr;
    memcpy(dest,left_ptr,left_ptr_max-left_ptr);
    return;
  }
  // General case
  while (left_ptr < left_ptr_max && right_ptr < right_ptr_max) {
    if (compar(right_ptr,left_ptr) < 0) {
      memcpy(dest,right_ptr,size);
      right_ptr+= size;
    }
    else {
      memcpy(dest,left_ptr,size);
      left_ptr+= size;
    }
    dest+= size;
  }
  if ( left_ptr < left_ptr_max )
    memcpy(dest,left_ptr,left_ptr_max-left_ptr);
  else if ( right_ptr < right_ptr_max )
    memcpy(dest,right_ptr,right_ptr_max-right_ptr);
  return;
}
size_t binary_search(void *base,void *val,size_t nmemb,size_t size,
		     int (*compar ) (const void *, const void * )){
  register size_t low_indx;
  register size_t high_indx;
  register size_t mid_indx;
  // special cases: val < array and val > array
  if ( compar(val,base) < 0 )
    return 0;
  else if ( compar(val,base+size*(nmemb-1) ) > 0 )
    return nmemb;
  low_indx= 0;
  high_indx= nmemb-1;
  while ( compar(base+low_indx*size,base+high_indx*size) < 0 ) {
    mid_indx= low_indx + (high_indx - low_indx) / 2;
    if ( compar(val,base+mid_indx*size) <= 0 )
      high_indx= mid_indx;
    else
      low_indx= mid_indx + 1;   
  }
  return low_indx;
}
void parallel_merge(void *src1, void *src2,
		    void *dest, 
		    size_t nmemb1, size_t nmemb2, size_t size, 
		    int (*compar ) (const void *, const void * )){
  register size_t val_indx;
  if ( nmemb1+nmemb2 < PARALLEL_SERIAL_MERGE_SWITCH ) {
    serial_merge(src1,src2,dest,nmemb1,nmemb2,size,compar);
    return;
  }
  if ( nmemb1 < nmemb2 ) {
    parallel_merge(src2,src1,dest,nmemb2,nmemb1,size,compar);
    return;
  }
  val_indx= binary_search(src2,src1+size*(nmemb1/2),nmemb2,size,compar);
  memcpy(dest+size*((nmemb1/2)+val_indx),src1+size*(nmemb1/2),size);
  #pragma omp task
  {
    parallel_merge(src1,src2,dest,nmemb1/2,val_indx,size,compar);
  }
  #pragma omp task
  {
    parallel_merge(src1+size*(nmemb1/2+1),src2+size*val_indx,
		   dest+size*((nmemb1/2+1)+val_indx),
		   nmemb1-nmemb1/2-1,nmemb2-val_indx,size,compar);
  }
  #pragma omp taskwait
}
void parallel_mergesort(void *src, void *dest,
			size_t nmemb, size_t size, 
			int (*compar ) (const void *, const void * ),
			bool src2dest){
  if ( nmemb < PARALLEL_SERIAL_SORT_SWITCH ) {
    if ( src2dest ) {
      qsort(src,nmemb,size,compar);
      memcpy(dest,src,nmemb*size);
    }
    else {
      qsort(dest,nmemb,size,compar);
      memcpy(src,dest,nmemb*size);
    }
    return;
  }
  #pragma omp task
  {
    parallel_mergesort(src,dest,nmemb/2,size,compar,!src2dest);
  }
  #pragma omp task
  {
    parallel_mergesort(src+size*(nmemb/2),dest+size*(nmemb/2),
		       nmemb-nmemb/2,size,compar,!src2dest);
  }
  #pragma omp taskwait
  if ( src2dest )
    parallel_merge(src,src+size*(nmemb/2),dest,nmemb/2,nmemb-nmemb/2,
		   size,compar);
  else
    parallel_merge(dest,dest+size*(nmemb/2),src,nmemb/2,nmemb-nmemb/2,
		   size,compar);
  return;
}
void parallel_sort(void *base, size_t nmemb, size_t size, 
		   int (*compar ) (const void *, const void * )){
  void* aux= malloc( nmemb * size );
  memcpy(aux,base,nmemb*size);
  #pragma omp parallel
  #pragma omp single
  parallel_mergesort(aux,base,nmemb,size,compar,true);
  free(aux);
  return;
}

