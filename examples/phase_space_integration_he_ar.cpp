#define USE_MPI
#include "hawaii.h"
#include "hep_hawaii.hpp"

#define HEAR_IMPLEMENTATION
#include "HeAr.h"

#define Rmin 6.597835932201344e+00

#define QP_SIZE 6

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    INIT_WRANK;
    INIT_WSIZE;

    uint32_t seed = 42; // mt_goodseed();

    pes = pes_lab;
    dipole = dipole_lab;

    double MU = m_He * m_Ar / (m_He + m_Ar); 
    MoleculeSystem ms = init_ms(MU, ATOM, ATOM, NULL, NULL, seed);

    CalcParams params = {};
    //params.ps                = PAIR_STATE_FREE_AND_METASTABLE;
    params.ps                = PAIR_STATE_ALL;
    params.sampler_Rmin      = 4.0;
    params.sampler_Rmax      = 50.0;
    params.initialM0_npoints = 10000000;
    params.pesmin            = V_HeAr(Rmin) / HTOCM;
    
    double T = 50.0;
    
    hep::mpi_vegas_callback<double>(hep::mpi_vegas_verbose_callback<double>);
    
    double pf_analytic = analytic_full_partition_function_by_V(&ms, T);
    PRINT0("analytic pf = %.5e\n", pf_analytic);

    double hep_ppf, hep_ppf_err;    
    mpi_perform_integration(&ms, integrand_pf, &params, T, 15, 10e6, &hep_ppf, &hep_ppf_err);
    double ppf_ratio = hep_ppf / pf_analytic;
    PRINT0("PPF ratio: %.5e\n\n", ppf_ratio);
   
    double hep_M0, hep_M0_err; 
    hep::mpi_vegas_callback<double>(hep::mpi_vegas_verbose_callback<double>);
    mpi_perform_integration(&ms, integrand_M0, &params, T, 15, 10e6, &hep_M0, &hep_M0_err);

    hep_M0     *= ZeroCoeff / pf_analytic;
    hep_M0_err *= ZeroCoeff / pf_analytic;
    PRINT0("HEP M0: %.5e +/- %.5e\n\n", hep_M0, hep_M0_err);

    // ------------------------------------------------------------------------------------
    params.partial_partition_function_ratio = ppf_ratio;

    double M0, M0_std;
    mpi_calculate_M0(&ms, &params, T, &M0, &M0_std); 

    PRINT0("M0 = %.10e +/- %.10e [%.10e ... %.10e]\n", M0, M0_std, M0-M0_std, M0+M0_std);
    PRINT0("Error: %.3f%%\n", M0_std/M0 * 100.0);

    PRINT0("\n\n\n");

    free_ms(&ms);
    
    MPI_Finalize();

    return 0;
}
