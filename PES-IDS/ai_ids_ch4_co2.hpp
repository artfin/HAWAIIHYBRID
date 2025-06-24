#ifndef AI_IDS_CH4CO2_H_
#define AI_IDS_CH4CO2_H_

#include <cmath>
#include <Eigen/Dense> 
#include <gsl/gsl_sf_legendre.h>

extern "C" {
    void dipole_init();
    void dipole_free();
    void dipole_vector(double *q, double dip[3]);
}

#endif // AI_IDS_CH4CO2_H_

