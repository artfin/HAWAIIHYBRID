#include "hawaii.h"
#include "angles_handler.hpp"
#include "ai_pes_co2ar_lib.hpp"
#include "ai_ids_co2ar_lib.hpp"

int main()
{
    uint32_t seed = 43; // mt_goodseed();
    
    pes_init();
    dipole_init();
    pes = pes_lab;
    dipole = dipole_lab;

    double MU = m_CO2 * m_Ar / (m_CO2 + m_Ar); 
    double I1[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE, ATOM, I1, NULL, seed);

    CalcParams params = {};
    params.ps                               = PAIR_STATE_FREE_AND_METASTABLE;
    params.sampler_Rmin                     = 4.5;
    params.sampler_Rmax                     = 40.0;
    params.initialM0_npoints                = 1000000;
    params.initialM2_npoints                = 2000000;
    params.partial_partition_function_ratio = 2.68854e+05;
    params.pesmin                           = -195.6337098547 / HTOCM;
    
    double T = 300.0;

    double M0, M0_std;
    calculate_M0(&ms, &params, T, &M0, &M0_std); 

    printf("M0 = %.10e +/- %.10e [%.10e ... %.10e]\n", M0, M0_std, M0-M0_std, M0+M0_std);
    printf("Error: %.3f%%\n", M0_std/M0 * 100.0);

    double M2, M2_std;
    calculate_M2(&ms, &params, T, &M2, &M2_std);
    
    printf("M2 = %.10e +/- %.10e [%.10e ... %.10e]\n", M2, M2_std, M2-M2_std, M2+M2_std);
    printf("Error: %.3f%%\n", M2_std/M2 * 100.0);

    // if (assert_float_is_equal_to(M0, 1.350e-09, 1e-11) > 0) {
    //     exit(1);
    // }


    free_ms(&ms);

    pes_free();

    return 0;
}
