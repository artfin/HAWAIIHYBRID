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

/*
 *  Copyright (C) 2026 A.Finenko & D.Chistikov 
 *  Distributed under the GNU General Public License, version 3
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */       

