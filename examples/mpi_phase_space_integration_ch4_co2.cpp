#define USE_MPI
#include "hawaii.h"
#include "hep_hawaii.hpp"

#include "ai_pes_ch4_co2_lib.hpp"
#include "ai_ids_ch4_co2_lib.hpp"

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    INIT_WRANK;
    INIT_WSIZE;

    uint32_t seed = 42; // mt_goodseed();
    
    pes_init();
    pes = pes_lab;
    dipole_init();
    dipole = dipole_lab;

    double MU = m_CH4 * m_CO2 / (m_CH4 + m_CO2); 
    double I1[3] = {II_CH4, II_CH4, II_CH4};
    double I2[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(MU, ROTOR, LINEAR_MOLECULE, I1, I2, seed);

    CalcParams params = {};
    params.ps                               = PAIR_STATE_BOUND;
    params.sampler_Rmin                     = 4.751;
    params.sampler_Rmax                     = 40.0;
    params.initialM0_npoints                = 100000;
    params.initialM2_npoints                = 100000;
    params.partial_partition_function_ratio = 1.0;
    params.pesmin                           = -342.934 / HTOCM;
    
    double T = 280.0;
    
    // ------------------------------------------------------------------------------------
    // -------------------------------  HEP  ----------------------------------------------
    // ------------------------------------------------------------------------------------
    hep::mpi_vegas_callback<double>(hep::mpi_vegas_verbose_callback<double>);
    
    double hep_M0, hep_M0_err; 
    mpi_perform_integration(&ms, integrand_M0, &params, T, 12, 1e2, &hep_M0, &hep_M0_err);

    double pf_analytic = analytic_full_partition_function_by_V(&ms, T);

    hep_M0     *= ZeroCoeff / pf_analytic;
    hep_M0_err *= ZeroCoeff / pf_analytic;
    PRINT0("HEP M0: %.5e +/- %.5e\n\n", hep_M0, hep_M0_err);

    double hep_M2, hep_M2_err;
    mpi_perform_integration(&ms, integrand_M2, &params, T, 9, 1e2, &hep_M2, &hep_M2_err);

    hep_M2     *= SecondCoeff / pf_analytic;
    hep_M2_err *= SecondCoeff / pf_analytic;
    PRINT0("HEP M2: %.5e +/- %.5e\n\n", hep_M2, hep_M2_err);

    double hep_ppf, hep_ppf_err;
    mpi_perform_integration(&ms, integrand_pf, &params, T, 15, 2e6, &hep_ppf, &hep_ppf_err);
    double ppf_ratio = hep_ppf / pf_analytic;
    PRINT0("PPF ratio: %.5e\n\n", ppf_ratio);
    assert(false); 
    // ------------------------------------------------------------------------------------
    
    params.partial_partition_function_ratio = ppf_ratio;

    double M0, M0_std;
    mpi_calculate_M0(&ms, &params, T, &M0, &M0_std); 

    PRINT0("M0 = %.10e +/- %.10e [%.10e ... %.10e]\n", M0, M0_std, M0-M0_std, M0+M0_std);
    PRINT0("Error: %.3f%%\n", M0_std/M0 * 100.0);
    
    if (assert_float_is_equal_to(M0, 8.29e-04, 1e-05) > 0) {
        MPI_Finalize();
        return 1; 
    }

    PRINT0("\n\n\n");

    double M2, M2_std;
    mpi_calculate_M2(&ms, &params, T, &M2, &M2_std); 
    PRINT0("M2 = %.10e +/- %.10e [%.10e ... %.10e]\n", M2, M2_std, M2-M2_std, M2+M2_std);
    PRINT0("Error: %.3f%%\n", M2_std/M2 * 100.0);
    
    if (assert_float_is_equal_to(M2, 5.37, 2e-2) > 0) {
        MPI_Finalize();
        return 1; 
    }


    free_ms(&ms);
    pes_free();

    MPI_Finalize();

    return 0;
}
