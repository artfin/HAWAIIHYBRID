/*
 * artfin 14.09.25
 *
 * This implementation is based on:
 *    Richard L. Burden & J. Douglas Faires, Numerical Analysis, Ninth Edition
 *    Algoithm 3.5 (page 173)
 *
 ************************************************************************
 * Reference Python program that uses Spline from Scipy.Interpolate
 *
 *  import numpy as np                                                 
 *  from scipy.interpolate import CubicSpline
 *  
 *  if __name__ == "__main__":
 *     x = np.array([-1.0, -0.5, 0.0, 0.5])
 *     y = np.array([0.86199480, 0.95802009, 1.0986123, 1.2943767])
 *  
 *     cs = CubicSpline(x, y, bc_type=((1, 0.155362), (1, 0.451863)))
 *  
 *     print("Coefficients:\n", cs.c)
 *
 *     x = -0.9
 *     print("x = {}, spline = {}".format(x, cs(x)))
 *
 *
 ************************************************************************
 * Example:
 *
 * #define CSPLINE_IMPLEMENTATION
 * #include "src/cspline.h"
 * 
 * int main()
 * {
 *     double x[] = {-1.0, -0.5, 0.0, 0.5};
 *     double y[] = {0.86199480, 0.95802009, 1.0986123, 1.2943767};
 *     double fp0 = 0.155362;
 *     double fpn = 0.451863;
 * 
 *     CSpline *spl = cspline_init(x, y, sizeof(x)/sizeof(x[0]), fp0, fpn);
 * 
 *     printf("i a[i] \t b[i] \t c[i] \t d[i]\n");
 *     for (size_t i = 0; i < spl->N; ++i) {
 *         printf("%zu \t %.6f \t %.6f \t %.6f \t %.6f\n", i, spl->a[i], spl->b[i], spl->c[i], spl->d[i]);
 *     }
 * 
 *     double yv;
 *     double xv = -0.9; 
 *     if (cspline_eval(spl, xv, &yv)) {
 *         printf("x = %.5e => y = %.10e\n", xv, yv);
 *     }
 * 
 *     cspline_free(spl);
 *     return 0;
 * }
 */

#ifndef CSPLINE_H_
#define CSPLINE_H_

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    size_t N;
    double *x;
    double *y;
    double fp0;
    double fpn;
    double *a;
    double *h;
    double *alpha;
    double *II;
    double *mu;
    double *z;
    double *b;
    double *c;
    double *d;
    double *mem;
} CSpline;

CSpline* cspline_init(double *x, double *y, size_t len, double fp0, double fpn); 
bool cspline_eval(CSpline *cs, double x, double *y);
void cspline_free(CSpline *cs);


#ifdef CSPLINE_IMPLEMENTATION

static const double CSPLINE_EPS = 1e-12;

CSpline* cspline_init(double *x, double *y, size_t len, double fp0, double fpn) 
{
    CSpline *cs = (CSpline*) malloc(1 * sizeof(CSpline)); 
    cs->N = len - 1;
    size_t n = cs->N;

    size_t capacity = 10*(n + 1) + n;
    double *mem = (double*) malloc(capacity * sizeof(double));
    assert((mem != NULL) && "ASSERT: cannot allocate enough memory\n");
    memset(mem, 0, capacity*sizeof(double));
    cs->mem = mem;
    
    cs->fp0 = fp0;
    cs->fpn = fpn;

    size_t offset = 0;
    cs->x = mem + offset;     offset += (n + 1);
    cs->y = mem + offset;     offset += (n + 1);
    cs->a = mem + offset;     offset += (n + 1);
    cs->h = mem + offset;     offset += n;
    cs->alpha = mem + offset; offset += (n + 1);
    cs->II = mem + offset;    offset += (n + 1);
    cs->mu = mem + offset;    offset += (n + 1);
    cs->z = mem + offset;     offset += (n + 1);
    cs->b = mem + offset;     offset += (n + 1);      
    cs->c = mem + offset;     offset += (n + 1);
    cs->d = mem + offset;     offset += (n + 1);

    memcpy(cs->x, x, (n + 1)*sizeof(double));
    memcpy(cs->y, y, (n + 1)*sizeof(double));
    memcpy(cs->a, y, (n + 1)*sizeof(double));
    
    for (size_t i = 0; i < n; ++i) {
        cs->h[i] = cs->x[i + 1] - cs->x[i];

        if (fabs(cs->h[i]) < CSPLINE_EPS) {
            printf("ERROR: x-values are too close to each other: difference x[%zu] and x[%zu] is less than expected epsilon\n",
                   i, i + 1);
            exit(1);
        }

        if (cs->h[i] < 0) {
            printf("ERROR: x-values have to be strictly increasing. Found x[%zu] = %.5e and x[%zu] = %.5e which must be increasing\n", 
                   i, cs->x[i], i + 1, cs->x[i + 1]);
            exit(1); 
        }
    }
    
    cs->alpha[0] = 3.0*(cs->a[1] - cs->a[0])/cs->h[0] - 3.0*fp0;
    cs->alpha[n] = 3.0*fpn - 3.0*(cs->a[n] - cs->a[n - 1])/cs->h[n - 1];

    for (size_t i = 1; i < n; ++i) {
        cs->alpha[i] = 3.0/cs->h[i]*(cs->a[i + 1]-cs->a[i]) - 3.0/cs->h[i - 1]*(cs->a[i]-cs->a[i - 1]);
    }
    
    cs->II[0] = 2.0 * cs->h[0];
    cs->mu[0] = 0.5;
    cs->z[0] = cs->alpha[0] / cs->II[0];

    for (size_t i = 1; i < n; ++i) {
        cs->II[i] = 2.0*(cs->x[i + 1] - cs->x[i - 1]) - cs->h[i - 1]*cs->mu[i - 1];
        cs->mu[i] = cs->h[i] / cs->II[i];
        cs->z[i] = (cs->alpha[i] - cs->h[i - 1]*cs->z[i - 1])/cs->II[i];
    }

    cs->II[n] = cs->h[n - 1]*(2.0 - cs->mu[n - 1]);
    cs->z[n] = (cs->alpha[n] - cs->h[n - 1]*cs->z[n - 1]) / cs->II[n];
    
    cs->c[n] = cs->z[n];

    for (int j = (int)n - 1; j >= 0; --j) {
        cs->c[j] = cs->z[j] - cs->mu[j] * cs->c[j + 1];
        cs->b[j] = (cs->a[j + 1] - cs->a[j])/cs->h[j] - cs->h[j]*(cs->c[j + 1] + 2.0*cs->c[j])/3.0;
        cs->d[j] = (cs->c[j + 1] - cs->c[j]) / 3.0 / cs->h[j];
    }

    return cs;
}

bool cspline_eval(CSpline *cs, double x, double *y)
{
    double xmin = cs->x[0];
    double xmax = cs->x[cs->N];
    if ((x < xmin) || (x > xmax)) {
        printf("ERROR: requested x = %.5e is outside interpolation range [%.5e -- %.5e]\n",
                x, xmin, xmax); 
        return false;
    }

    size_t j = 0;
    while ((j <= cs->N - 1) && (x >= cs->x[j + 1])) j++; 

    *y = cs->a[j] + cs->b[j]*(x - cs->x[j]) + cs->c[j]*(x - cs->x[j])*(x - cs->x[j]) + cs->d[j]*(x - cs->x[j])*(x - cs->x[j])*(x - cs->x[j]); 
    return true; 
}

void cspline_free(CSpline *cs) 
{
    free(cs->mem);
    free(cs);
}

#endif // CSPLINE_IMPLEMENTATION

#endif // CSPLINE_H_
