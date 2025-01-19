#ifndef AI_PES_CO2AR_H_
#define AI_PES_CO2AR_H_

#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_sf_legendre.h>

#include "../constants.h"

#ifdef __cplusplus
extern "C" {
#endif

void init_pes();
void free_pes();

double pes_co2ar(double R, double Theta);
void dpes_co2ar(double R, double Theta, double *dR, double *dTheta);

#ifdef __cplusplus
}
#endif

// minimum parameters
// double RMIN = 6.4966051519;
// double THETAMIN = M_PI / 2.0;
// double PESMIN = -195.6337098547; // cm-1 

#endif // AI_PES_CO2AR_H_
