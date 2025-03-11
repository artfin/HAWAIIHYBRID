#define USE_MPI
#include "hawaii.h"
#include "hep_hawaii.hpp"

#define HEAR_IMPLEMENTATION
#include "HeAr.h"

#define Rmin 6.597835932201344e+00

double pes(double *q) {
    return V_HeAr(q[2]);
} 

void dpes(double *q, double *dpesdq) {
    dpesdq[0] = 0.0;           // dVdPhi
    dpesdq[1] = 0.0;           // dVdTheta
    dpesdq[2] = dV_HeAr(q[2]); // dVdR
}

void dipole_lab(double *q, double diplab[3]) {
    double Phi    = q[0];
    double Theta  = q[1];
    double diplen = dip_HeAr(q[2]);

    diplab[0] = diplen * sin(Theta) * cos(Phi); 
    diplab[1] = diplen * sin(Theta) * sin(Phi);
    diplab[2] = diplen * cos(Theta);
}

#define QP_SIZE 6

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    INIT_WRANK;
    INIT_WSIZE;

    uint32_t seed = 42; // mt_goodseed();

    dipole = dipole_lab;

    double MU = m_He * m_Ar / (m_He + m_Ar); 
    MoleculeSystem ms = init_ms(MU, ATOM, ATOM, NULL, NULL, seed);

    CalcParams params = {};
    params.ps           = FREE_AND_METASTABLE;
    params.sampler_Rmin = 4.0;
    params.sampler_Rmax = 40.0;
    params.pesmin       = V_HeAr(Rmin) / HTOCM;
    
    double T = 300.0;
    
    double pf_analytic = analytic_full_partition_function_by_V(&ms, T);
    
    hep::mpi_vegas_callback<double>(hep::mpi_vegas_verbose_callback<double>);

    double hep_ppf, hep_ppf_err;    
    mpi_perform_integration(&ms, integrand_pf, &params, T, 15, 2e6, &hep_ppf, &hep_ppf_err);
    double ppf_ratio = hep_ppf / pf_analytic;
    PRINT0("PPF ratio: %.5e\n\n", ppf_ratio);

    free_ms(&ms);
    
    MPI_Finalize();

    return 0;
}
