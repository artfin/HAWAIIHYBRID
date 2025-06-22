#ifndef AI_PES_CO2AR_HPP_
#define AI_PES_CO2AR_HPP_

extern "C" {
    void pes_init();
    void pes_free();
    double pes_lab(double *q);
    void dpes_lab(double *q, double *dpesdq);
}

#endif // AI_PES_CO2AR_HPP_

