#ifndef LOESS_H_
#define LOESS_H_

#include <iostream>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include <omp.h>
#include <float.h>

/*
 * LOESS: locally weighted polynomial regression
 *
 * LOESS blends the simplicity of linear least squares regression
 * with the adaptability of nonlinear regression. It achieves this by fitting simple
 * model to localized subset of the data, gradually constructing a function that captures 
 * the deterministic pattern of the variation in the data -- effectively filtering out 
 * the random component that follows some probability distribution.
 *
 * Degree of local polynomials:
 * The local polynomials fitted to each subset of the data are typically of either first
 * or second degree. Employing a zero-degree polynomial reduces LOESS to a weighted moving 
 * average. While higher-degree polynomials could theoretically be used, they are not aligned 
 * with the spirit of LOESS. Such polynomials are prone to overfitting within each subset and 
 * thus often lead to numerical instability. 
 *
 * Weight function:
 * The weight function, gives the moset weight to the data points nearest the point
 * of estimation and the least weight to the data points that are furthest away.
 * The use of the weights is based on the idea that points near each other are more
 * likely to be related to each other in a simple way than points that are further
 * apart.
 * The traditional weight function used for LOESS is the tricube weight function:
 * w(x) = (1 - |x|^3)^3 for |x| < 1 and 0 otherwise.
 * The main criteria for the weight function are the following (Cleveland, 1979):
 *  - w(x) > 0 for |x| < 1 since negative weights do not make sense
 *  - w(-x) = w(x): there is no reason to treat points to the left of x differently
 *    from those to the right 
 *  - w(x) is a nonincreasing function for x >= 0: it seems unreasonable to allow
 *    a point that is closer to x to have less weight than the one that is further
 *    away 
 *  - w(x) = 0 for |x| >= 1
 * In addition it seems desirable that w(x) decrease smoothly to 0 as x goes from 0
 * to 1. Such a weight function is more likely to produce a smoothed result. The
 * tricube has been chosen since it enhances a chi-squared distributional approximation
 * of an estimate of the error variance. So it should provide an adequate smooth in 
 * many situations.
 * The weight for a specific point in any localized subset of data is obtained by 
 * evaluating the weight function at the distance between that point and the point 
 * of estimation, after scaling the distance so that the maximum absolute distance 
 * over all of the points in the subset of data is exactly one.
 *
 * Blog post with python implementation:
 * https://towardsdatascience.com/loess-373d43b03564
 * The method is proposed by William S. Cleveland in 1979. 
 * https://www.tandfonline.com/doi/abs/10.1080/01621459.1979.10481038  
 *
 * Apparently, LOESS is equivalent to Savitzky-Golay filtering.
 *
 * NIST example of LOESS Computation:  
 * https://www.itl.nist.gov/div898/handbook/pmd/section1/dep/dep144.htm 
 */

typedef struct {
    size_t *items;
    size_t count;
    size_t capacity;
} Window;

extern bool loess_debug;

typedef enum {
    WEIGHT_TRICUBE,
    WEIGHT_BISQUARE,
} WEIGHT_FUNC;

typedef enum {
    LS_COMPLETE_ORTHOGONAL_DECOMPOSITION,
    LS_CHOLESKY_SOLVER,
    LS_QR, 
} LS_METHOD; 

typedef struct {
    size_t degree;   // degree of local polynomial [recommended: 2-3]
    size_t ws_min;   // minimum window size 
    double ws_step;  // window size step
    size_t ws_delay; // optional: the index at which the window is starting to increase
    size_t ws_cap;   // optional: cap on window size 
} Smoothing_Config;

extern WEIGHT_FUNC loess_weight;

void loess_init(double *x, double *y, size_t len);
double loess_estimate(double x, size_t window_size, size_t degree);
double *loess_create_grid(double xmin, double xmax, size_t npoints);
double *loess_apply_smoothing(Smoothing_Config *config);
// double *apply_smoothing(size_t degree, double xmin, double xmax, size_t npoints, size_t wsmin, double wsstep, size_t wsdelay);
void loess_free();

#endif // LOESS_H_
