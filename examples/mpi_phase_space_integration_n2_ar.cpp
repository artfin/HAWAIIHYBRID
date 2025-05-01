#define USE_MPI
#include "hawaii.h"
#include "hep_hawaii.hpp"

#include "angles_handler.hpp"

#include "ai_pes_n2_ar_pip_nn.hpp"
#include "ai_ids_n2_ar_pip_nn.hpp"

double pes(double *q)
// PHI THETA R PHI1T THETA1T  
{
    XYZ xyz = lab_to_xyz(q);
    return pes_xyz(&xyz); 
} 

void dpes(double *q, double *dq) {
    UNUSED(q);
    UNUSED(dq);
    TODO("dpes");
}

void dipole_lab(double *q, double diplab[3]) {
    double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);
    
    Eigen::Vector3d dipmol = dipole_bf(qmol[0], qmol[4]); 
    
    double sinphiem, cosphiem;
    double sinthetaem, costhetaem;
    double sinpsiem, cospsiem;

    sincos(qmol[1], &sinphiem, &cosphiem);
    sincos(qmol[2], &sinthetaem, &costhetaem);
    sincos(qmol[3], &sinpsiem, &cospsiem);

    Sz_filler(Sphiem, sinphiem, cosphiem);
    Sx_filler(Sthetaem, sinthetaem, costhetaem);
    Sz_filler(Spsiem, sinpsiem, cospsiem);
       
    Eigen::Vector3d diplab_eig = Sphiem.transpose() * Sthetaem.transpose() * Spsiem.transpose() * dipmol; 
   
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
    
    pes_init(false); 
    dipole_init(false);
    dipole = dipole_lab;

    double MU = m_N2 * m_Ar / (m_N2 + m_Ar); 
    double I1[2] = {II_N2, II_N2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE, ATOM, I1, NULL, seed);

    CalcParams params = {};
    params.ps                               = FREE_AND_METASTABLE;
    params.sampler_Rmin                     = 4.5;
    params.sampler_Rmax                     = 40.0;
    params.initialM0_npoints                = 10000000;
    params.initialM2_npoints                = 10000000;
    params.partial_partition_function_ratio = 1.0;
    params.pesmin                           = -97.3 / HTOCM;
    
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
    
    //if (assert_float_is_equal_to(M0, 3.64e-04, 8e-6) > 0) {
    //    MPI_Finalize();
    //    return 1; 
    //}

    PRINT0("\n\n\n");

    double M2, M2_std;
    mpi_calculate_M2(&ms, &params, T, &M2, &M2_std); 
    PRINT0("M2 = %.10e +/- %.10e [%.10e ... %.10e]\n", M2, M2_std, M2-M2_std, M2+M2_std);
    PRINT0("Error: %.3f%%\n", M2_std/M2 * 100.0);
    
    //if (assert_float_is_equal_to(M2, 6.405e-01, 2e-2) > 0) {
    //    MPI_Finalize();
    //    return 1; 
    //}

    free_ms(&ms);

    MPI_Finalize();

    return 0;
}
