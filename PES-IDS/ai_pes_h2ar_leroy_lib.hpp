#ifndef AR_H2_PES_LEROY_DERIVATIVES_LIB_H_
#define AR_H2_PES_LEROY_DERIVATIVES_LIB_H_

extern "C" {
    double pes_lab(double *q);
    void dpes_lab(double *q, double *dpesdq);
}

#endif // AR_H2_PES_LEROY_DERIVATIVES_LIB_H_

