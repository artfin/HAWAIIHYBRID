#include "array.h"

void init_array(Array *a, double *data, size_t n) {
    assert(a->data != NULL);
    assert(a->n >= n);
    memcpy(a->data, data, n * sizeof(double));
}

Array create_array(size_t n) {
    Array a = {0};
    a.data = (double*) malloc(n * sizeof(double));
    assert(a.data != NULL);
    a.n = n; 
    return a;
}

Array arena_create_array(Arena *a, size_t n) {
    Array arr = {0};
    arr.data = (double*) arena_alloc(a, n * sizeof(double));
    assert(arr.data != NULL);
    arr.n = n;
    return arr;
}

void free_array(Array *a) {
    free(a->data);
}

void print_array(Array a) {
    printf("Array[%zu] = ", a.n); 
    for (size_t i = 0; i < a.n; ++i) {
        printf("%.3f ", a.data[i]);
    }
    printf("\n"); 
}


