#define USE_MPI
#include "hawaii.h"

#include "array.h"
#include "trajectory.h"

#include "ai_pes_ch4_co2_lib.hpp"
#include "ai_ids_ch4_co2_lib.hpp"

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    INIT_WRANK;
    INIT_WSIZE;

    int seeds[] = {20106943, 74537365, 68838871, 45228894, 87924815, 75409565, 78079463, 71551420, 61270358, 12907113, 37812429, 43846516, 78563783, 67020198, 46463763, 15913615, 56260494, 73487105, 39711609, 17700181, 26290757, 10277061, 79863361, 29916933, 49969708, 56772087, 50637984, 10326549, 38641292, 81352997, 39721216, 73377986, 25410770, 65662596, 95180457, 12378249, 90265129, 87941783, 28243843, 40892493, 60538578, 67701455, 54511007, 95983369, 82062787, 59498255, 74798224, 29737064, 62225324, 18939487, 99787768, 87346783, 85009677, 63086756, 14433074, 99516596, 75192796, 75650862, 35503295, 20831228, 33162524, 12437412, 91210209, 76571171, 17922249, 15052896, 20751049, 49130285, 47619032, 29471832, 37775141, 74794610, 29101489, 55424490, 18716635, 20356657, 32172602, 97293267, 77861181, 29092863, 62639840, 56987692, 81112748, 58656852, 27469169, 27154689, 49220900, 95017234, 62779012, 97324386, 96580635, 45610890, 41955768, 14529933, 46207757, 44406420, 65851255, 34626188, 40063013, 66832671};
    assert(_wsize <= (int) (sizeof(seeds)/sizeof(int)) && "ASSERT: more seeds should be added\n");

    double MU = m_CH4 * m_CO2 / (m_CH4 + m_CO2); 
    double I1[3] = {II_CH4, II_CH4, II_CH4};
    double I2[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(MU, ROTOR, LINEAR_MOLECULE, I1, I2, seeds[_wrank]);
    
    pes_init();
    pes = pes_lab;
    dpes = dpes_lab;

    dipole_init();
    dipole = dipole_lab;

    double tolerance = 1e-12;
    Trajectory traj = init_trajectory(&ms, tolerance);
   
    CalcParams params = {};
    params.ps                               = PAIR_STATE_FREE_AND_METASTABLE;
    params.sampler_Rmin                     = 4.751;
    params.sampler_Rmax                     = 40.0;
    params.niterations                      = 3;
    params.total_trajectories               = 600;
    params.cvode_tolerance                  = 1e-12;
    params.sampling_time                    = 200.0;
    params.MaxTrajectoryLength              = 65536;
    params.Rcut                             = 40.0;
    params.partial_partition_function_ratio = 2.69441e+05;
    params.initialM0_npoints                = 2000;
    params.initialM2_npoints                = 2000;
    params.pesmin                           = -342.934 / HTOCM;
    params.cf_filename                      = "./CF-CH4-CO2-F-300.0.txt";

    printf("MU = %.15e\n", MU);
    printf("II(CH4) = %.15e\n", II_CH4);
    printf("pesmin = %.15e\n", params.pesmin);

    assert(false);

    double Temperature = 300.0;
    
    CFnc cf = calculate_correlation_and_save(&ms, &params, Temperature);

    if (_wrank == 0) {
        printf("\n\n");
        printf("Correlation function is calculated. The initial values are:\n");
        for (size_t i = 0; i < 10; ++i) {
            printf("%.2f   %.10e\n", cf.t[i], cf.data[i]*ALU*ALU*ALU);
        }
    }

    free_trajectory(&traj);
    free_ms(&ms);
    free_cfnc(cf); 

    MPI_Finalize();

    return 0;
}
