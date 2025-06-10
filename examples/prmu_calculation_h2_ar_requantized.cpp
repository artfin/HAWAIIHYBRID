#define JINI_HISTOGRAM_BINS 11
#define JINI_HISTOGRAM_MAX  5.0 
#define JFIN_HISTOGRAM_BINS 11
#define JFIN_HISTOGRAM_MAX  5.0 
#define NSWITCH_HISTOGRAM_BINS 10
#define NSWITCH_HISTOGRAM_MAX 10.0
#define USE_MPI

#include "hawaii.h"

#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

#include "ai_pes_h2ar_leroy.h"
#include "ai_ids_h2_ar_pip_nn.hpp"

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    Arena a = {};
    
    double MU = m_H2 * m_Ar / (m_H2 + m_Ar); 
    double I1[2] = {II_H2, II_H2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE_REQ_HALFINTEGER, ATOM, I1, NULL, 0);

    printf("II_H2: %.15e\n", II_H2);
    printf("MU: %.15e\n", MU);
    printf("D_H2: %.15e\n", D_H2);

    Array qp0 = arena_create_array(&a, ms.QP_SIZE);
    get_qp_from_ms(&ms, &qp0);

    for (size_t i = 0; i < ms.QP_SIZE; ++i) {
        printf("%.5e\n", qp0.data[i]);
    }
    
    double l = l_H2*sqrt(5.6919435017e+01 / 60.8530119); 
    double II = (m_H/2.0*l*l); 
    double B_cm = Planck/(8.0*M_PI*M_PI*II*AMU*ALU*ALU) / LightSpeed_cm; // cm-1
    printf("II_H2 = %.10e => II = %.10e\n", II_H2, II);
    printf("B_cm = %.10e\n", B_cm);
    printf("l = %.10e\n", l);

    MPI_Finalize();

    return 0;
}

int main2(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    INIT_WRANK;
    INIT_WSIZE;

    int seeds[] = {20106943, 74537365, 68838871, 45228894, 87924815, 75409565, 78079463, 71551420, 61270358, 12907113, 37812429, 43846516, 78563783, 67020198, 46463763, 15913615, 56260494, 73487105, 39711609, 17700181, 26290757, 10277061, 79863361, 29916933, 49969708, 56772087, 50637984, 10326549, 38641292, 81352997, 39721216, 73377986, 25410770, 65662596, 95180457, 12378249, 90265129, 87941783, 28243843, 40892493, 60538578, 67701455, 54511007, 95983369, 82062787, 59498255, 74798224, 29737064, 62225324, 18939487, 99787768, 87346783, 85009677, 63086756, 14433074, 99516596, 75192796, 75650862, 35503295, 20831228, 33162524, 12437412, 91210209, 76571171, 17922249, 15052896, 20751049, 49130285, 47619032, 29471832, 37775141, 74794610, 29101489, 55424490, 18716635, 20356657, 32172602, 97293267, 77861181, 29092863, 62639840, 56987692, 81112748, 58656852, 27469169, 27154689, 49220900, 95017234, 62779012, 97324386, 96580635, 45610890, 41955768, 14529933, 46207757, 44406420, 65851255, 34626188, 40063013, 66832671};
    assert(_wsize <= (int) (sizeof(seeds)/sizeof(int)) && "ASSERT: more seeds should be added\n");
    
    double MU = m_H2 * m_Ar / (m_H2 + m_Ar); 
    double I1[2] = {II_H2, II_H2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE_REQ_HALFINTEGER, ATOM, I1, NULL, seeds[_wrank]);
    ms.m1.DJ = D_H2;
    
    dipole_init(true);
    dipole = dipole_lab;

    CalcParams params = {};
    params.ps = FREE_AND_METASTABLE;
    // TODO:
    // used in q_generator to generate phase-point for pr/mu and then overwritten, which should probably be corrected
    // however, these values are used to obtain phase-space estimates for M0/M2, so they are needed anyway  
    params.sampler_Rmin = 4.0;
    params.sampler_Rmax = 40.0;

    params.niterations                      = 10;
    params.total_trajectories               = 100;
    params.cvode_tolerance                  = 1e-8;
    params.sampling_time                    = 200.0;
    params.MaxTrajectoryLength              = 65536;
    params.initialM0_npoints                = 10000;
    params.initialM2_npoints                = 10000;
    params.pesmin                           = -105.0 / HTOCM;
    params.R0                               = 40.0;
    params.ApproximateFrequencyMax          = 1500.0;
    params.torque_cache_len                 = 30;
    params.torque_limit                     = 5e-6;
    params.jini_histogram_bins              = 101;
    params.jini_histogram_max               = 10.0;
    params.jfin_histogram_bins              = 101;
    params.jfin_histogram_max               = 10.0;
    params.odd_j_spin_weight                = 0.75;
    params.even_j_spin_weight               = 0.25;
    
    double Temperature = 300.0;

    String_Builder sf_filename = {};
    sb_append_format(&sf_filename, "./SF-PRMU-H2-Ar-%.1f.txt", Temperature); 
    params.sf_filename = sf_filename.items; 

    SFnc sf = calculate_spectral_function_using_prmu_representation_and_save(&ms, &params, Temperature); 
    
    if (_wrank == 0) {
        printf("\n\n");
        printf("Spectral function is calculated. The initial values are:\n");
        for (size_t i = 0; i < 10; ++i) {
            printf("%.2f   %.10e\n", sf.nu[i], sf.data[i]);
        }
    }

    sb_free(&sf_filename);
    free_ms(&ms);
    free_sfnc(sf); 
     
    MPI_Finalize();

    return 0;
}

