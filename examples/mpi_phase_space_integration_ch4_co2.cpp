#define USE_MPI
#include "hawaii.h"
#include "hep_hawaii.hpp"

#include "ai_pes_ch4_co2.h"
#include "ai_ids_ch4_co2.hpp"

#include "angles_handler.hpp"

#define NKAL 5

double pes(double *qlab) {
    static double qkal[NKAL];
    CH4_linear_molecule_lab_to_kal(qlab, qkal); 
    return pes_ch4co2(qkal[0], qkal[1], qkal[2], qkal[3], qkal[4]); 
} 

void dpes(double *q, double *dq) {
    UNUSED(q);
    UNUSED(dq);
    TODO("dpes");
}

void dipole_lab(double *qlab, double diplab[3]) {
    double qkal[NKAL];
    CH4_linear_molecule_lab_to_kal(qlab, qkal);
    
    double dipkal[3];
    dipole_vector(qkal, dipkal);

    double sinphi, cosphi;
    double sintheta, costheta;
    double sinpsi, cospsi;

    sincos(qlab[3], &sinphi, &cosphi);
    sincos(qlab[4], &sintheta, &costheta);
    sincos(qlab[5], &sinpsi, &cospsi);

    Sz_filler(Sphi1t, sinphi, cosphi);
    Sx_filler(Stheta1t, sintheta, costheta);
    Sz_filler(Spsi1t, sinpsi, cospsi);
      
    Eigen::Vector3d dipkal_eig = Eigen::Map<Eigen::Vector3d>(dipkal, 3);
    Eigen::Vector3d diplab_eig = Sphi1t.transpose() * Stheta1t.transpose() * Spsi1t.transpose() * dipkal_eig; 
   
    diplab[0] = diplab_eig(0); 
    diplab[1] = diplab_eig(1); 
    diplab[2] = diplab_eig(2); 
}



int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    INIT_WRANK;
    INIT_WSIZE;

    uint32_t seed = 42; // mt_goodseed();
    
    init_pes();
    init_ids();
    dipole = dipole_lab;

    double MU = m_CH4 * m_CO2 / (m_CH4 + m_CO2); 
    double I1[3] = {II_CH4, II_CH4, II_CH4};
    double I2[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(MU, ROTOR, LINEAR_MOLECULE, I1, I2, seed);

    CalcParams params = {};
    params.ps                               = FREE_AND_METASTABLE;
    params.sampler_Rmin                     = 4.751;
    params.sampler_Rmax                     = 40.0;
    params.initialM0_npoints                = 10000000;
    params.initialM2_npoints                = 10000000;
    params.partial_partition_function_ratio = 1.0;
    params.pesmin                           = -342.934 / HTOCM;
    
    double T = 300.0;
    
    // ------------------------------------------------------------------------------------
    // -------------------------------  HEP  ----------------------------------------------
    // ------------------------------------------------------------------------------------
    double hep_M0, hep_M0_err; 
    hep::mpi_vegas_callback<double>(hep::mpi_vegas_verbose_callback<double>);
    mpi_perform_integration(&ms, integrand_M0, &params, T, 12, 1e6, &hep_M0, &hep_M0_err);

    double pf_analytic = analytic_full_partition_function_by_V(&ms, T);

    hep_M0     *= ZeroCoeff / pf_analytic;
    hep_M0_err *= ZeroCoeff / pf_analytic;
    PRINT0("HEP M0: %.5e +/- %.5e\n\n", hep_M0, hep_M0_err);

    double hep_M2, hep_M2_err;
    mpi_perform_integration(&ms, integrand_M2, &params, T, 9, 1e6, &hep_M2, &hep_M2_err);

    hep_M2     *= SecondCoeff / pf_analytic;
    hep_M2_err *= SecondCoeff / pf_analytic;
    PRINT0("HEP M2: %.5e +/- %.5e\n\n", hep_M2, hep_M2_err);

    double hep_ppf, hep_ppf_err;
    mpi_perform_integration(&ms, integrand_pf, &params, T, 12, 1e6, &hep_ppf, &hep_ppf_err);
    double ppf_ratio = hep_ppf / pf_analytic;
    PRINT0("PPF ratio: %.5e\n\n", ppf_ratio); 
    // ------------------------------------------------------------------------------------
    
    params.partial_partition_function_ratio = ppf_ratio;

    double M0, M0_std;
    mpi_calculate_M0(&ms, &params, T, &M0, &M0_std); 

    PRINT0("M0 = %.10e +/- %.10e [%.10e ... %.10e]\n", M0, M0_std, M0-M0_std, M0+M0_std);
    PRINT0("Error: %.3f%%\n", M0_std/M0 * 100.0);
    
    if (assert_float_is_equal_to(M0, 3.64e-04, 8e-6) > 0) {
        MPI_Finalize();
        return 1; 
    }

    PRINT0("\n\n\n");

    double M2, M2_std;
    mpi_calculate_M2(&ms, &params, T, &M2, &M2_std); 
    PRINT0("M2 = %.10e +/- %.10e [%.10e ... %.10e]\n", M2, M2_std, M2-M2_std, M2+M2_std);
    PRINT0("Error: %.3f%%\n", M2_std/M2 * 100.0);
    
    if (assert_float_is_equal_to(M2, 6.405e-01, 2e-2) > 0) {
        MPI_Finalize();
        return 1; 
    }


    free_ms(&ms);
    free_pes();

    MPI_Finalize();

    return 0;
}
