#include "loess.hpp"

#include "hawaii.h"

static double *xp = NULL;
static double *yp = NULL;
static size_t len = 0;
static double minx = FLT_MAX;
static double maxx = FLT_MIN;
static double miny = FLT_MAX; 
static double maxy = FLT_MIN;

bool loess_debug = false;
WEIGHT_FUNC loess_weight = WEIGHT_TRICUBE;

#undef da_append
#undef da_insert

// @TODO: introduce more mechanisms of fiddling with 'eval_with_linear_window_size'
// at least we should cap the window size at some value
//
// @TODO: test loess with valgrind, we definitely leak some memory in 'evaluate'
//
// @TODO: document individual functions and copy the docstrings to pdf
//
// @TODO: OpenMP parallelization of 'eval_with_linear_window_size' function
//
// @TODO: test different weight functions 
//
// we use C++ compiler for this translation unit and it has some peculiar problem
// of not being able to implicitly cast void* to the type of (da)->items which is 
// size_t* in this file. I decided to hardcode the cast in the macro... 
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


void init_loess(double *x, double *y, size_t ilen) 
{
    xp = (double*) malloc(ilen * sizeof(double));
    yp = (double*) malloc(ilen * sizeof(double));
    len = ilen;

    for (size_t i = 0; i < len; ++i) {
        if (x[i] > maxx) maxx = x[i];
        if (x[i] < minx) minx = x[i];
        if (y[i] > maxy) maxy = y[i];
        if (y[i] < miny) miny = y[i];
    }

    for (size_t i = 0; i < len; ++i) {
        xp[i] = (x[i] - minx) / (maxx - minx);
        yp[i] = (y[i] - miny) / (maxy - miny);
    }
} 

Window make_window(double* distances, size_t window_size) 
{
    size_t min_idx = 0;
    double min_val = distances[0];

    for (size_t i = 1; i < len; ++i) {
        double val = distances[i];
        if (val < min_val) {
            min_val = val;
            min_idx = i;
        }
    }

    if (min_idx == 0) {
        return (Window) {
            .items = linspace_size_t(0, window_size - 1, window_size),
            .count = window_size,
            .capacity = window_size,
        };
    } else if (min_idx == len - 1) {
        return (Window) {
            .items= linspace_size_t(len - window_size, len - 1, window_size),
            .count = window_size,
            .capacity = window_size,
        };
    }

    Window window{};
    da_append(&window, min_idx);
    
    while (window.count < window_size) {
        size_t ind_lower = window.items[0]; // current lower boundary of the window
        size_t ind_upper = window.items[window.count - 1]; // current upper boundary of the window

        if (ind_lower == 0) {
            da_append(&window, ind_upper + 1);
        } else if (ind_upper == len - 1) {
            da_insert(&window, 0, ind_lower - 1);
        } else if (distances[ind_lower - 1] < distances[ind_upper + 1]) {
            da_insert(&window, 0, ind_lower - 1);
        } else {
            da_append(&window, ind_upper + 1);
        }
    }
    
    return window;
}

inline double tricube(double x) {
    return (1.0 - x * x * x) * (1.0 - x * x * x) * (1.0 - x * x * x);  
}

inline double bisquare(double x) {
    return (1.0 - x * x) * (1.0 - x * x);  
}

double estimate(double x, size_t window_size, size_t degree)
{
    double nx = (x - minx) / (maxx - minx);
    
    double* distances = (double*) malloc(len * sizeof(double));
    for (size_t i = 0; i < len; ++i) {
        distances[i] = fabs(xp[i] - nx);
    }
    
    Window window = make_window(distances, window_size);
    if (loess_debug) {
        printf("DEBUG: window: ");
        for (size_t i = 0; i < window.count; ++i) {
            printf("%zu ", window.items[i]);
        }
        printf("\n");
    }


    double max_distance = FLT_MIN;
    for (size_t i = 0; i < window.count; ++i) {
        if (distances[window.items[i]] > max_distance) max_distance = distances[window.items[i]]; 
    }
    
    /* calculate weights */    
    Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(window.count, window.count); 
    
    if (loess_debug) printf("DEBUG: wm:\n");
    for (size_t i = 0; i < window.count; ++i) {
        double window_distance = distances[window.items[i]];
        switch (loess_weight) {
            case WEIGHT_TRICUBE: {
              weights(i, i) = tricube(window_distance / max_distance);
              break;
            }
            case WEIGHT_BISQUARE: {
              weights(i, i) = bisquare(window_distance / max_distance);
              break;
            }
        }
        if (loess_debug) printf("%.4e ", weights(i, i));
    }

    if (loess_debug) printf("\n");

   
    // @naming
    Eigen::MatrixXd xm = Eigen::MatrixXd::Zero(window.count, degree+1);

    if (loess_debug) printf("DEBUG: xm:\n");

    for (size_t i = 0; i < window.count; ++i) {
        double xv = xp[window.items[i]];

        for (size_t j = 0; j < degree+1; ++j) {
            xm(i, j) = pow(xv, j);
            if (loess_debug) printf("%.4e ", xm(i, j));
        }
        if (loess_debug) printf("\n");
    }

    Eigen::MatrixXd xmt_wm = xm.transpose() * weights;
    
    if (loess_debug) {
        printf("DEBUG: xmt_wm:\n");
        for (size_t i = 0; i < degree+1; ++i) {
            for (size_t j = 0; j < window.count; ++j) {
                printf("%.4e ", xmt_wm(i, j));
            }
            printf("\n");
        }
    }
    
    Eigen::VectorXd yv = Eigen::VectorXd::Zero(window.count, 1);
    for (size_t i = 0; i < window.count; ++i) {
        yv(i) = yp[window.items[i]];
    }
 
    // note: this is hard to reproduce in GSL so I refused the idea to move this code completely to C
    // and, more importantly, Eigen is more optimised so it does not make much sense to switching back
    // to GSL/BLAS  
    Eigen::VectorXd beta = (xmt_wm * xm).completeOrthogonalDecomposition().pseudoInverse() * xmt_wm * yv;
    if (loess_debug) {
        printf("DEBUG: beta: ");
        for (size_t i = 0; i < degree+1; ++i) {
            printf("%.4e ", beta(i));
        }
        printf("\n");
    }

    Eigen::VectorXd x_pow = Eigen::VectorXd::Zero(degree + 1);
    for (size_t i = 0; i < degree + 1; ++i) {
        x_pow(i) = std::pow(nx, i);
    }

    double yn = beta.dot(x_pow);
    double y = yn * (maxy - miny) + miny; 
    
    if (loess_debug) { 
        printf("DEBUG: estimate: x = %.5e, y = %.5e, window_size = %zu, degree = %zu\n", x, y, window_size, degree); 
    }

    return y;
}

double *eval_with_linear_window_size(size_t degree, double xmin, double xmax, size_t npoints, size_t wsmin, double wsstep, size_t wsdelay, double xraw_step)
{
    double *smoothed = (double*) malloc(npoints * sizeof(double));
    memset(smoothed, 0.0, npoints * sizeof(double));

    double xstep = 0;
    if (npoints > 1) {
        xstep = (xmax - xmin) / (npoints - 1);
    }

    size_t window_size = wsmin;

    for (size_t i = 0; i < npoints; ++i) {
        double x = xmin + xstep * i;
       
        if (i > wsdelay) window_size = (size_t)(wsmin + wsstep*i);

        if (window_size < 1) {
            printf("ERROR: at iteration %zu the window size does not contain any points. Exiting...\n", i);
            return smoothed;
        } 

        smoothed[i] = estimate(x, window_size, degree);

        if (i % 100 == 0) {
            printf("INFO: i = %zu, frequency = %.3e, smoothed value = %.3e, window_size (in points) = %zu, window_size (in frequency units) = %.3e\n",
                    i, x, smoothed[i], window_size, window_size * xraw_step); 
        }
    } 

    return smoothed;
}

