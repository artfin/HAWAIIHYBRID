#ifndef AI_IDS_CH4CO2_H_
#define AI_IDS_CH4CO2_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_sf_legendre.h>

#include <Eigen/Dense>

void init_ids();
void free_ids();
void dipole_vector(double *q, double dip[3]);

#endif // AI_IDS_CH4CO2_H_

