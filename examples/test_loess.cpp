#include <stdio.h>

#include "loess.hpp"

#define N 6

#define DA_INIT_CAP 256

#define da_append(da, item)                                                                    \
    do {                                                                                       \
        if ((da)->count >= (da)->capacity) {                                                   \
            (da)->capacity = (da)->capacity == 0 ? DA_INIT_CAP : (da)->capacity*2;             \
            (da)->items = (size_t*) realloc((da)->items, (da)->capacity*sizeof(*(da)->items)); \
            assert((da)->items != NULL && "ASSERT: not enough memory\n");                      \
        }                                                                                      \
                                                                                               \
        (da)->items[(da)->count++] = (item);                                                   \
    } while (0)

#define da_insert(da, i, item)                                                                       \
    do {                                                                                             \
        if ((i < 0) || ((i) > (da)->count)) {                                                        \
            assert(0 && "ASSERT: index out of bounds\n");                                            \
        }                                                                                            \
        if ((da)->count >= (da)->capacity) {                                                         \
            (da)->capacity = (da)->capacity == 0 ? DA_INIT_CAP : (da)->capacity*2;                   \
            (da)->items = (size_t*) realloc((da)->items, (da)->capacity*sizeof(*(da)->items));       \
            assert((da)->items != NULL && "ASSERT: not enough memory\n");                            \
        }                                                                                            \
        memmove((da)->items + (i) + 1, (da)->items + (i), ((da)->count - (i))*sizeof(*(da)->items)); \
        (da)->items[(i)] = (item);                                                                   \
        (da)->count++;                                                                               \
    } while(0)

#define da_last(da) (da)->items[(da)->count - 1]

extern "C" {
  double pes(double *q) { (void) q; return 0; }
  void dpes(double *q, double *dVdq) { (void)q; (void)dVdq; }
}

void test_da_array() {
    Window window{};
    da_append(&window, 1);
    da_append(&window, 2);
    da_append(&window, 3);

    da_insert(&window, 0, 4);

    printf("window: ");
    for (size_t i = 0; i < window.count; ++i) {
        printf("%zu ", window.items[i]);
    }
    printf("\n");

    printf("last: %zu\n", da_last(&window));
}

int main()
{
    test_da_array();

    double x[N] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double y[N] = {0.1, 0.2, 0.3, 0.4, 0.6, 0.9}; 

    printf("-- NEW VERSION\n");
    loess_init(x, y, N);
    double u = loess_estimate(2.5, 3 /*window size*/, 3);
    printf("u = %.10e\n", u);

    return 0;
}
