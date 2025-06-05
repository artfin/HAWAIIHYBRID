#define USE_MPI
#include "hawaii.h"

#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

double pes(double *q) {
    double r = q[0];

    if (r < 6.0) { 
        return 100.0;
    }

	return 0.0;
} 

void dpes(double *q, double *dpesdq) {
    (void) q;
    memset(dpesdq, 0, 10*sizeof(double));    
}

void dipole_lab(double *q, double diplab[3]) 
// PHI THETA R PHI1T THETA1T  
{
    diplab[0] = mu_CO * sin(q[4]) * cos(q[3]); 
    diplab[1] = mu_CO * sin(q[4]) * sin(q[3]);
    diplab[2] = mu_CO * cos(q[4]);

#if 0
    double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);
    
    double dipmol[3];
    dipmol[0] = arco_dipx_ind(qmol[0] /* R */, qmol[4] /* Theta */);
    dipmol[1] = 0.0; 
    dipmol[2] = arco_dipz_ind(qmol[0] /* R */, qmol[4] /* Theta */);

    double sinphiem, cosphiem;
    double sinthetaem, costhetaem;
    double sinpsiem, cospsiem;

    sincos(qmol[1], &sinphiem, &cosphiem);
    sincos(qmol[2], &sinthetaem, &costhetaem);
    sincos(qmol[3], &sinpsiem, &cospsiem);

    Sz_filler(Sphiem, sinphiem, cosphiem);
    Sx_filler(Sthetaem, sinthetaem, costhetaem);
    Sz_filler(Spsiem, sinpsiem, cospsiem);
       
    Eigen::Vector3d dipmol_eig = Eigen::Map<Eigen::Vector3d>(dipmol, 3);
    Eigen::Vector3d diplab_eig = Sphiem.transpose() * Sthetaem.transpose() * Spsiem.transpose() * dipmol_eig; 
   
    diplab[0] = diplab_eig(0); 
    diplab[1] = diplab_eig(1); 
    diplab[2] = diplab_eig(2);
#endif
}


int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    INIT_WRANK;
    INIT_WSIZE;

    int seeds[] = {20106943, 74537365, 68838871, 45228894, 87924815, 75409565, 78079463, 71551420, 61270358, 12907113, 37812429, 43846516, 78563783, 67020198, 46463763, 15913615, 56260494, 73487105, 39711609, 17700181, 26290757, 10277061, 79863361, 29916933, 49969708, 56772087, 50637984, 10326549, 38641292, 81352997, 39721216, 73377986, 25410770, 65662596, 95180457, 12378249, 90265129, 87941783, 28243843, 40892493, 60538578, 67701455, 54511007, 95983369, 82062787, 59498255, 74798224, 29737064, 62225324, 18939487, 99787768, 87346783, 85009677, 63086756, 14433074, 99516596, 75192796, 75650862, 35503295, 20831228, 33162524, 12437412, 91210209, 76571171, 17922249, 15052896, 20751049, 49130285, 47619032, 29471832, 37775141, 74794610, 29101489, 55424490, 18716635, 20356657, 32172602, 97293267, 77861181, 29092863, 62639840, 56987692, 81112748, 58656852, 27469169, 27154689, 49220900, 95017234, 62779012, 97324386, 96580635, 45610890, 41955768, 14529933, 46207757, 44406420, 65851255, 34626188, 40063013, 66832671};
    assert(_wsize <= (int) (sizeof(seeds)/sizeof(int)) && "ASSERT: more seeds should be added\n");

    double MU = m_CO * m_Ar / (m_CO + m_Ar); 
    double I1[2] = {II_CO, II_CO};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE_REQ_INTEGER, ATOM, I1, NULL, seeds[_wrank]);
    ms.m1.DJ = D_CO;

    dipole = dipole_lab;

    CalcParams params = {};
    params.ps = FREE_AND_METASTABLE;
    // TODO:
    // used in q_generator to generate phase-point for pr/mu and then overwritten, which should probably be corrected
    // however, these values are used to obtain phase-space estimates for M0/M2, so they are needed anyway  
    params.sampler_Rmin = 4.0;
    params.sampler_Rmax = 40.0;

    params.niterations                                   = 10;
    params.total_trajectories                            = 100;
    params.cvode_tolerance                               = 1e-8;
    params.sampling_time                                 = 200.0;
    params.MaxTrajectoryLength                           = 2097152;
    params.allow_truncating_trajectories_at_length_limit = false;
    params.partial_partition_function_ratio              = 2.68854e+05;
    params.initialM0_npoints                             = 10000;
    params.initialM2_npoints                             = 10000;
    params.pesmin                                        = -105.0 / HTOCM;
    params.R0                                            = 40.0;
    params.ApproximateFrequencyMax                       = 150.0;
    params.jini_histogram_bins                           = 351;
    params.jini_histogram_max                            = 35.0;
    params.jfin_histogram_bins                           = 351;
    params.jfin_histogram_max                            = 35.0;
    params.torque_cache_len                              = 30;
    params.torque_limit                                  = 5e-6;
    params.average_time_between_collisions               = 0.988e-10 / ATU; // a.t.u.
    
    double Temperature = 300.0;
   
    String_Builder sf_filename = {};
    sb_append_format(&sf_filename, "./SF-PRMU-line-test-%.1f.txt", Temperature); 
    params.sf_filename = sf_filename.items; 

    SFnc sf = calculate_spectral_function_using_prmu_representation_and_save(&ms, &params, Temperature); 
    
    if (_wrank == 0) {
        printf("\n\n");
        printf("Spectral function is calculated. The initial values are:\n");
        for (size_t i = 0; i < 10; ++i) {
            printf("%.4f   %.10e\n", sf.nu[i], sf.data[i]);
        }
    }

    sb_free(&sf_filename);
    free_ms(&ms);
    free_sfnc(sf); 
     
    MPI_Finalize();

    return 0;
}

