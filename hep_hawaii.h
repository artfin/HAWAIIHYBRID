#ifndef HAWAII_HEP_H_
#define HAWAII_HEP_H_

#include "hawaii.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    INTEGRAND_M0,
    INTEGRAND_M2,
    INTEGRAND_PF,
} IntegrandType;

void c_mpi_perform_integration(MoleculeSystem *ms, IntegrandType integrand_type, CalcParams *params, double Temperature, size_t niterations, size_t npoints, double *m, double *q);

#ifdef __cplusplus
}
#endif

#endif // HAWAII_HEP_H_
