#ifndef DYNAMICAL_ALLOCATIONS_H
#define DYNAMICAL_ALLOCATIONS_H

#include <stdlib.h>
#include <stdio.h>

inline void my_malloc(size_t** ptr,size_t s) 
{
  *ptr = malloc(sizeof(size_t) + s);
  *ptr[0] = s;
}

inline void my_realloc(size_t ** ptr, size_t s)
{
  *ptr = realloc(*ptr, sizeof(size_t) + s);
  *ptr[0] = s;
}

inline void my_free(void * ptr) 
{
  free(ptr);
}

inline size_t allocated_size(void * ptr) 
{
  return ((size_t*)ptr)[0];
}

inline size_t size_t_array_length(void * ptr) 
{
  return ((size_t*)ptr)[0]/sizeof(size_t);

}

#endif
