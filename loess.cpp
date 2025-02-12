#include "loess.hpp"

#include "hawaii.h"

static double *xp = NULL;
static double *yp = NULL;
static size_t len = 0;
static double minx = FLT_MAX;
static double maxx = FLT_MIN;
static double miny = FLT_MAX; 
static double maxy = FLT_MIN;

static double XRAW_STEP = 0;

static double GRID_XMIN = FLT_MAX;
static double GRID_XMAX = FLT_MIN;
static size_t GRID_NPOINTS = 0;
static double GRID_XSTEP = 0;

bool loess_debug = false;
WEIGHT_FUNC loess_weight = WEIGHT_TRICUBE;
LS_METHOD ls_method = LS_COMPLETE_ORTHOGONAL_DECOMPOSITION;

#undef da_append
#undef da_insert

// @TODO: check if 'thread_local' is needed anywhere
// @TODO: document individual functions and copy the docstrings to pdf
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

#define da_last(da) (da)->items[(da)->count - 1]


#ifndef _OPENMP
int omp_get_num_threads() { return 1; }
int omp_get_thread_num() { return 1; }
#endif // _OPENMP

void loess_init(double *x, double *y, size_t ilen)
/*
 * This function prepares the input data for LOESS by storing the predictor 'x' and response 'y' values,
 * and setting up the necessary values for subsequent LOESS calculations. The function assumes that
 * the input arrays `x` and `y` are of equal length and contain valid numerical data. 
 */
{
    xp = (double*) malloc(ilen * sizeof(double));
    yp = (double*) malloc(ilen * sizeof(double));
    len = ilen;

    if (len < 2) {
        printf("ERROR: expected length of arrays to be at least 2\n");
        exit(1);
    }

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
    
    XRAW_STEP = x[1] - x[0];
} 

void loess_free()
{
    free(xp);
    free(yp);
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


    size_t ind_lower = min_idx;
    size_t ind_upper = (min_idx == len - 1) ? min_idx : min_idx + 1;
    for (size_t window_count = 0; window_count < window_size; window_count++) {
        if (ind_lower == 0) {
            ind_upper++;
        } else if (ind_upper == len - 1) {
            ind_lower--;
        } else if (distances[ind_lower - 1] < distances[ind_upper + 1]) {
            ind_lower--;
        } else {
            ind_upper++;
        }
    }
    
    Window window{
        .items    = (size_t*) malloc(window_size * sizeof(double)),
        .count    = window_size,
        .capacity = window_size,
    };

    for (size_t i = 0; i < window_size; ++i) {
        window.items[i] = ind_lower + i; 
    }

    // @note:
    // this seems to be an extremely dynamic way of accumulating the window
    // maybe we can rewrite this thing to avoid considering elements one-by-one 
    // and thus getting rid of constant appending/inserting?
    /* 
    Window window{};
    da_append(&window, min_idx);
  
    while (window.count < window_size) {
        size_t ind_lower = window.items[0]; // current lower boundary of the window
        size_t ind_upper = da_last(&window); // current upper boundary of the window
        printf("ind_lower: %zu, ind_upper: %zu\n", ind_lower, ind_upper);

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
    */    
        
    return window;
}

inline double tricube(double x) {
    return (1.0 - x * x * x) * (1.0 - x * x * x) * (1.0 - x * x * x);  
}

inline double bisquare(double x) {
    return (1.0 - x * x) * (1.0 - x * x);  
}

double loess_estimate(double x, size_t window_size, size_t degree)
{
    double nx = (x - minx) / (maxx - minx);
    
    double* distances = (double*) malloc(len * sizeof(double));
    for (size_t i = 0; i < len; ++i) {
        distances[i] = fabs(xp[i] - nx);
    }

    // @note: 
    // the indices of the window are expected to be in the increasing order    
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

    // the weighted linear regression problem is formulated as
    // (X^T W X) beta = X^T W y
    
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(window.count, degree+1);

    if (loess_debug) printf("DEBUG: xm:\n");

    for (size_t i = 0; i < window.count; ++i) {
        double xv = xp[window.items[i]];

        for (size_t j = 0; j < degree+1; ++j) {
            X(i, j) = pow(xv, j);
            if (loess_debug) printf("%.4e ", X(i, j));
        }
        if (loess_debug) printf("\n");
    }

    Eigen::MatrixXd XtW = X.transpose() * weights;
    
    if (loess_debug) {
        printf("DEBUG: xmt_wm:\n");
        for (size_t i = 0; i < degree+1; ++i) {
            for (size_t j = 0; j < window.count; ++j) {
                printf("%.4e ", XtW(i, j));
            }
            printf("\n");
        }
    }
    
    Eigen::VectorXd Y = Eigen::VectorXd::Zero(window.count, 1);
    for (size_t i = 0; i < window.count; ++i) {
        Y(i) = yp[window.items[i]];
    }

    // discussion of approaches to solve linear squares systems in Eigen:
    // https://eigen.tuxfamily.org/dox/group__LeastSquares.html 
    //
    // beta = (X^T W X)^(-1) X^T W Y 
    // 
    // @note: this is hard to reproduce in GSL so I refused the idea to move this code completely to C
    // and, more importantly, Eigen is more optimised so it does not make much sense to switching back
    // to GSL/BLAS  
    Eigen::VectorXd beta;
    switch (ls_method) {
        case LS_COMPLETE_ORTHOGONAL_DECOMPOSITION: {
            // Option 1: Complete orthogonal decomposition 
            // @timing: CH4-CO2, 300K, 11.40m -> 10.30m 
            beta = (XtW * X).completeOrthogonalDecomposition().pseudoInverse() * XtW * Y;
            break; 
        }
        case LS_CHOLESKY_SOLVER: {
            // Option 2: normal equations
            // this thing produces more or less the same result with approximately the same speed
            beta = (XtW * X).ldlt().solve(XtW * Y); 
            break;
        }
        case LS_QR_NO_PIVOTING: {
            // Option 3: QR decomposition // 11.51m
            beta = (XtW * X).householderQr().solve(XtW * Y);
            break;
        } 
    } 

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
    
    free(distances);
    free(window.items);
    
    return y;
}

double *loess_create_grid(double grid_xmin, double grid_xmax, size_t grid_npoints) 
{
    if (grid_xmin < minx) {
        printf("ERROR: the starting point of the grid (%.5e) is less than the minimum xvalue in the raw data (%.5e). Adjust the grid bounds to be within the data range.\n", 
                grid_xmin, minx);
        return NULL;
    }

    if (grid_xmax > maxx) {
        printf("ERROR: the ending point of the grid (%.5e) is greater than the maximum xvalue in the raw data (%.5e). Adjust the grid bounds to be within the data range.\n",
                grid_xmax, maxx);
        return NULL;
    }

    if (grid_xmin > grid_xmax) {
        printf("ERROR: invalid grid range: the starting point (%.5e) is greater than the ending point (%.5e). The grid must be defined such that grid_xmin < grid_xmax.\n",
                grid_xmin, grid_xmax);
        return NULL;
    }

    if (grid_npoints < 1) {
        printf("ERROR: invalid number of grid points: %zu. The grid must contain > 1 points\n", grid_npoints);
        return NULL;
    }
    
    GRID_XMIN = grid_xmin;
    GRID_XMAX = grid_xmax;
    GRID_NPOINTS = grid_npoints; 

    GRID_XSTEP = (GRID_XMAX - GRID_XMIN) / (GRID_NPOINTS - 1);

    return linspace(GRID_XMIN, GRID_XMAX, GRID_NPOINTS);
}

double *loess_apply_smoothing(Smoothing_Config *config) 
{
    if (GRID_NPOINTS <= 0) {
        printf("ERROR: the grid is not initialized. Ensure 'loess_create_grid' is called to define the grid before applying the smoothing algorithm\n");
        return NULL;
    }

    double *smoothed = (double*) malloc(GRID_NPOINTS * sizeof(double));
    memset(smoothed, 0.0, GRID_NPOINTS * sizeof(double));

    bool should_exit = false;

    #pragma omp parallel
    {
        printf("INFO: loess_apply_smoothing is run using %d threads\n", omp_get_num_threads());
    }

    #pragma omp parallel for schedule(dynamic, 50)
    for (size_t i = 0; i < GRID_NPOINTS; ++i) {
        if (should_exit) continue;
        double x = GRID_XMIN + GRID_XSTEP * i;
      
        size_t window_size = config->ws_min; 
        if (i > config->ws_delay) window_size = (size_t)(config->ws_min + config->ws_step*i);

        if (config->ws_cap > 0) {
            if (window_size > config->ws_cap) window_size = config->ws_cap; 
        }

        if (window_size < 1) {
            printf("ERROR: at iteration %zu the window size does not contain any points. Exiting...\n", i);
            should_exit = true; 
        } 

        smoothed[i] = loess_estimate(x, window_size, config->degree);

        if (i % 100 == 0) {
            printf("(%d) INFO: i = %zu, frequency = %.3e, smoothed value = %.3e, window_size (in points) = %zu, window_size (in frequency units) = %.3e\n",
                    omp_get_thread_num(), i, x, smoothed[i], window_size, window_size * XRAW_STEP); 
        }
    } 

    return smoothed;
}
