
#ifndef AI_PES_CH4CO2_H_
#define AI_PES_CH4CO2_H_

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <gsl/gsl_sf_legendre.h>

#include "../constants.h"

#ifdef __cplusplus
extern "C" {
#endif

void pes_init();
void pes_free();

double pes_ch4co2(double R, double phi1, double theta1, double phi2, double theta2);
void dpes_ch4co2(double R, double phi1, double theta1, double phi2, double theta2, double *dR, double *dphi1, double *dtheta1, double *dphi2, double *dtheta2);

double dpesdR(double R, double phi1, double theta1, double phi2, double theta2);
double dpesdphi1(double R, double phi1, double theta1, double phi2, double theta2);
double dpesdtheta1(double R, double phi1, double theta1, double phi2, double theta2);
double dpesdphi2(double R, double phi1, double theta1, double phi2, double theta2);
double dpesdtheta2(double R, double phi1, double theta1, double phi2, double theta2);

inline size_t gsl_index(int l, int m) {
    return l * (l + 1) / 2 + m;
}

#ifdef __cplusplus
}
#endif 

// minimum parameters
// double RMIN = 6.4966051519;
// double THETAMIN = M_PI / 2.0;
// double PESMIN = -195.6337098547; // cm-1 

#endif // AI_PES_CH4CO2_H_
