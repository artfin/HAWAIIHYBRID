#define USE_MPI
#include "hawaii.h"

#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

// В работе Pedersen написано, что потенциал предполагает, что theta_m = 0 отвечает C-O-Ar, а theta_m = Pi  —  O-C-Ar, где theta_m определен как угол между вектором от C к О и вектором, 
// соединяющим центр масс CO и Ar. В лабораторной системе координат мы вправе двумя способами выбрать углы (phi1t, theta1t). В ситуации theta1t = 0 молекула ориентирована вдоль оси OZ и 
// положительную координату имеет или  атом C или атом O: (+)C — O(-) или (-) С — O (+). Пусть мы считаем, что theta1t = 0 отвечает ситуации (-) C — O (+). 
// То есть phi1t = (произвольное значение), theta1t = 0, Phi = (произвольное значение), Theta = 0 обозначает конфигурацию Ar-O-C, где все атомы лежат на оси OZ лабораторной системы координат. 
// Применяя пересчет углов в подвижную систему, получаем theta_m = 0, что согласуется с конвенцией, принятой в потенциале Pederson'a.

// Dipole Rizzo
// θ = 0 deg corresponds to Ar-O-C an θ = 180 deg corresponds Ar-C-O
// Dipole direction: -CO+

extern "C" {
	// IPOT = 3 : Use the parameters of the CCSD(T)/aug-cc-pVQZ-33211
	//            3D surface at R1. This surface is probably not
	//            trustworthy beyond the interval (this is NOT tested!)
	//               1.898 bohr <= R1 <= 2.234 bohr
	//            due to the single-reference nature of CCSD(T).
	//            Thus, intended only for 2-D calculations at
	//            a given CO distance in the above trust interval.

    // R1    - CO distance in bohr
    // R2    - distance from CO center-of-mass to Ar in bohr
    // XCOS2 - cos to the intermolecular angle
	void potv(double* res, double* r1, double* r2, double* xcos2);
	void potv_d(double* v, double* vd, double* r1, double* r1d, double* r2, double* r2d,
		        double* xcos2, double* xcos2d);
}

#define DIPOLE_COAR_IMPLEMENTATION
#include "dipole_coar.c"

double pes(double *q) {
    static double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);

    // r1 - CO bond length in bohr
    // r2 - distance from CO center-of-mass to Ar in bohr
    // theta - angle in rad
	double _r1 = l_CO;
	double _r2 = qmol[0];
	double _xcos2 = cos(qmol[4]);

    double r;
	potv(&r, &_r1, &_r2, &_xcos2);
	return r;
} 

void dpes(double *q, double *dpesdq) {
    static Eigen::Matrix<double, 5, 5> jac;
    static Eigen::Matrix<double, 5, 1> derivatives_mol, derivatives_lab; 
    static double qmol[5];
     
    jac.setZero();
    linear_molecule_atom_lab_to_mol(q, qmol);
    linear_molecule_atom_Jacobi_mol_by_lab(jac, q, qmol);  
    
	double _r1 = l_CO; 
	double _r1d = 0.0;
	double _r2 = qmol[0]; 
	double _r2d = 1.0;
	double _xcos2 = cos(qmol[4]);
	double _xcos2d = 0.0;

	double t, dR, dTheta;
	
    potv_d(&t, &dR, &_r1, &_r1d, &_r2, &_r2d, &_xcos2, &_xcos2d);

	_r1d = 0.0;
	_r2d = 0.0;
	_xcos2d = 1.0;
	potv_d(&t, &dTheta, &_r1, &_r1d, &_r2, &_r2d, &_xcos2, &_xcos2d);

    // [dpes_dR, 0, 0, 0, dpes_dTheta]
    derivatives_mol(0) = dR; 
    derivatives_mol(4) = -sin(qmol[4])*dTheta; 

    //for (size_t i = 0; i < 5; ++i) {
    //    printf("derivatives_mol(%zu) = %.10lf\n", i, derivatives_mol(i));
    //}
    //
    //for (size_t i = 0; i < 5; ++i) {
    //    for (size_t j = 0; j < 5; ++j) {
    //        printf("%.5e ", jac(i, j)); 
    //    }
    //    printf("\n");
    //}
    
    derivatives_lab = jac * derivatives_mol;
    Eigen::VectorXd::Map(dpesdq, 5) = derivatives_lab;
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

int main2()
{
    double B_MHz = Planck/(8.0*M_PI*M_PI*II_CO*AMU*ALU*ALU) / 1e6; // MHz 
    double B_cm = Planck/(8.0*M_PI*M_PI*II_CO*AMU*ALU*ALU) / LightSpeed_cm; // cm-1
    printf("B(CO) = %.5e cm-1\n", B_cm);
    printf("B(CO) = %.5f MHz\n", B_MHz);

    double mu_D = mu_CO * ADIPTODEBYE;
    printf("mu(CO) = %.5f D\n", mu_D); 

    return 0;
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
    params.MaxTrajectoryLength              = 2097152; 
    params.partial_partition_function_ratio = 2.68854e+05;
    params.initialM0_npoints                = 10000;
    params.initialM2_npoints                = 10000;
    params.pesmin                           = -105.0 / HTOCM;
    params.R0                               = 40.0;
    params.ApproximateFrequencyMax          = 150.0;
    params.torque_cache_len                 = 30;
    params.torque_limit                     = 5e-6;
    params.average_time_between_collisions  = 0.988e-10 / ATU; // a.t.u.
   
    String_Builder sf_filename = {};
    sb_append_format(&sf_filename, "./SF-PRMU-CO-Ar-300.0-POISSON.txt"); 
    params.sf_filename = sf_filename.items; 

    double Temperature = 300.0;
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

