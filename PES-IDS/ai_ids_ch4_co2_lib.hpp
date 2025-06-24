#ifndef AI_IDS_CH4_CO2_LIB_HPP_
#define AI_IDS_CH4_CO2_LIB_HPP_

#include <Eigen/Dense>

extern "C" {
    void dipole_init();
    void dipole_free();
    void dipole_lab(double *qlab, double diplab[3]);
}


#endif // AI_IDS_CH4_CO2_LIB_HPP_

