#ifndef ARRAY_H_
#define  ARRAY_H_

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

typedef struct {
    double *data;
    size_t n;
} Array;

#ifdef __cplusplus
extern "C" {
#endif

Array create_array(size_t n);
void free_array(Array *a);

void init_array(Array *a, double *data, size_t n);

void print_array(Array a);

#ifdef __cplusplus
}
#endif

#endif //  ARRAY_H_

