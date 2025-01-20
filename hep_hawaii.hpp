#ifndef HAWAII_HEP_
#define HAWAII_HEP_

#include "hawaii.h"

// included in <complex.h> and collides with templates in hep
#undef I

#include <hep/mc-mpi.hpp>

typedef double (*Integrand)(hep::mc_point<double> const&);

void transform_variables(hep::mc_point<double> const& x, double* transformed, double* Jac);

double integrand_M0(hep::mc_point<double> const& x);
double integrand_M2(hep::mc_point<double> const& x);
double integrand_pf(hep::mc_point<double> const& x);

void mpi_perform_integration(MoleculeSystem *ms, Integrand integrand, CalcParams *params, double T, size_t niterations, size_t npoints, double *m, double *q);


#endif // HAWAII_HEP_
