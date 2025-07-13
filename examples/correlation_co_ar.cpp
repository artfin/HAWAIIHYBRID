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
#include "dipole_coar.cpp"

double pes_lab(double *q) {
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
    //printf("r = %.3e, theta = %.3e => V = %.3e cm-1\n", _r2, qmol[4], r*HTOCM);
	return r;
} 

void dpes_lab(double *q, double *dpesdq) {
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

const char *var_to_cstring(int n) {
    switch (n) {
        case 0: return "Phi";
        case 1: return "pPhi";
        case 2: return "Theta";
        case 3: return "pTheta";
        case 4: return "R";
        case 5: return "pR";
        case 6: return "Phi1";
        case 7: return "pPhi1";
        case 8: return "Theta1";
        case 9: return "pTheta1";
    }

    UNREACHABLE("var_to_cstring");
}

void test_rhs(MoleculeSystem *ms, Array qp)
{
    printf("\n-----------------------------------------\n");
    printf("Testing analytic derivatives of Hamiltonian against numerical ones\n");
    printf("The derivatives are shown in the same order used in 'rhs' function.\n");

    put_qp_into_ms(ms, qp);

    N_Vector y    = make_vector(ms->QP_SIZE);
    N_Vector ydot = make_vector(ms->QP_SIZE);

    memcpy(N_VGetArrayPointer(y), qp.data, qp.n * sizeof(double));
    rhs(0.0, y, ydot, (void*) ms);    

    size_t order = 6;
    Array num_derivatives = compute_numerical_rhs(ms, order);

    printf("# \t analytic \t numeric \t difference \n");
    for (size_t i = 0; i < ms->QP_SIZE; ++i) {
        printf("dot(%s): %.10e \t %.10e \t %.10e\n", var_to_cstring(i), NV_Ith_S(ydot, i), num_derivatives.data[i], NV_Ith_S(ydot, i) - num_derivatives.data[i]);

        if (assert_float_is_equal_to(NV_Ith_S(ydot, i), num_derivatives.data[i], 1e-4) > 0) {
            exit(1);
        }
    }
    printf("-----------------------------------------\n\n\n");
    
    put_qp_into_ms(ms, qp);
    
    free_array(&num_derivatives); 
    N_VDestroy(y);
    N_VDestroy(ydot);
}


void test_trajectory()
{
    int seed = 42;

    double MU = m_CO * m_Ar / (m_CO + m_Ar); 
    double I1[2] = {II_CO, II_CO};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE, ATOM, I1, NULL, seed);
    
    pes = pes_lab;
    dpes = dpes_lab;
    dipole = dipole_lab;

    Array qp = create_array(ms.QP_SIZE);
    double data[] = {
        1.748408280024098e-01,
        8.000086892719947e+00,
        2.043243374007463e-01,
        1.292311668112794e+01,
        2.966821669009366e+01,
       -3.786729452372973e-00,
        1.462013919064840e+00,
        1.199991310728005e+01,
        1.700437238862704e+00,
        3.160851788770303e+01,
    }; 
    
    init_array(&qp, data, ms.QP_SIZE);
    put_qp_into_ms(&ms, qp);

    test_rhs(&ms, qp);

    double tolerance = 1e-12;
    Trajectory traj = init_trajectory(&ms, tolerance);
   
    double E0 = Hamiltonian(&ms);
    set_initial_condition(&traj, qp);
    print_array(qp);

    CalcParams params = {};
    params.sampling_time = 200.0; 
    
    double t = 0.0;
    double tout = params.sampling_time;

    double dipt[3];

    for (size_t nstep = 0; ; ++nstep, tout += params.sampling_time)
    {
        int status = make_step(&traj, tout, &t);

        if (status > 0) {
            fprintf(stderr, "CVODE ERROR: status = %d\n", status);
            exit(1);
        }
      
        // 21.12.2024 NOTE: 
        // We copy the "N_Vector y" from "Trajectory" into "MoleculeSystem" on each call of rhs.
        // However after "rhs" has been called, CVode makes a step on dynamic variables. So we have 
        // to update the phase-point in MoleculeSystem after "make_step" function has returned. 
        put_qp_into_ms(&ms, (Array){.data = N_VGetArrayPointer(traj.y), .n = ms.QP_SIZE});
        
        double E = Hamiltonian(&ms);

        extract_q_and_write_into_ms(&ms);
        dipole_lab(ms.intermediate_q, dipt);

        printf("%10.1lf \t %12.10lf \t %10.12lf \t %.10f %.10f %.10f\n", t, ms.intermolecular_qp[IR], E-E0, dipt[0], dipt[1], dipt[2]);

        if (ms.intermolecular_qp[IR] > 30.0) break;

        //if (assert_float_is_equal_to(E-E0, 0.0, 1e-8) > 0) {
        //    exit(1); 
        //}
    }

    free_trajectory(&traj);
    free_ms(&ms);
    free_array(&qp);
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
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE, ATOM, I1, NULL, seeds[_wrank]);

    dipole_1 = dipole_lab;
    dipole_2 = dipole_CO_lab;
    pes = pes_lab;
    dpes = dpes_lab;
    
    double tolerance = 1e-12;
    Trajectory traj = init_trajectory(&ms, tolerance);
   
    CalcParams params = {};
    params.ps                               = PAIR_STATE_FREE_AND_METASTABLE;
    params.sampler_Rmin                     = 4.0;
    params.sampler_Rmax                     = 30.0;
    params.niterations                      = 3;
    params.total_trajectories               = 600;
    params.cvode_tolerance                  = 1e-12;
    params.sampling_time                    = 200.0;
    params.MaxTrajectoryLength              = 65536;
    params.Rcut                             = 40.0;
    params.partial_partition_function_ratio = 2.68854e+05;
    params.initialM0_npoints                = 1000000;
    params.initialM2_npoints                = 1000000;
    params.pesmin                           = -105.0 / HTOCM;
    params.cf_filename                      = "./CF-CO-Ar-F-300.0.txt";
    
    printf("MU = %.5e\n", MU);
    printf("II_CO = %.5e\n", II_CO);
    printf("pesmin = %.5e\n", params.pesmin);


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
    free(ms.m1.torque_cache);

    MPI_Finalize();

    return 0;
}

