#ifndef AI_IDS_CO2AR_LIB_HPP_
#define AI_IDS_CO2AR_LIB_HPP_

extern "C" {
    void dipole_init();
    void dipole_lab(double *q, double diplab[3]);
}

#endif // AI_IDS_CO2AR_LIB_HPP_
