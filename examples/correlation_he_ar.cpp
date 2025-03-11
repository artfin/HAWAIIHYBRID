#define USE_MPI
#include "hawaii.h"

#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

#define HEAR_IMPLEMENTATION
#include "HeAr.h"

double pes(double *q) {
    return V_HeAr(q[2]);
} 

void dpes(double *q, double *dpesdq) {
    dpesdq[0] = 0.0;           // dVdPhi
    dpesdq[1] = 0.0;           // dVdTheta
    dpesdq[2] = dV_HeAr(q[2]); // dVdR
}

void dipole_lab(double *q, double diplab[3]) {
    double Phi    = q[0];
    double Theta  = q[1];
    double diplen = dip_HeAr(q[2]);

    diplab[0] = diplen * sin(Theta) * cos(Phi); 
    diplab[1] = diplen * sin(Theta) * sin(Phi);
    diplab[2] = diplen * cos(Theta);
}

const char *var_to_cstring(int n) {
    switch (n) {
        case 0: return "Phi";
        case 1: return "pPhi";
        case 2: return "Theta";
        case 3: return "pTheta";
        case 4: return "R";
        case 5: return "pR";
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

    size_t order = 2;
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

    double MU = m_He * m_Ar / (m_He + m_Ar); 
    MoleculeSystem ms = init_ms(MU, ATOM, ATOM, NULL, NULL, seed);

    Array qp = create_array(ms.QP_SIZE);
    double data[] = {
        1.748408280024098e-01,
        8.000086892719947e+00,
        2.043243374007463e-01,
        1.292311668112794e+01,
        2.966821669009366e+01,
       -3.786729452372973e-00,
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
      
        double E = Hamiltonian(&ms);

        extract_q_and_write_into_ms(&ms);
        dipole_lab(ms.intermediate_q, dipt);

        printf("%10.1lf \t %12.10lf \t %12.15lf %.10f %.10f %.10f\n", t, ms.intermolecular_qp[IR], E-E0, dipt[0], dipt[1], dipt[2]);

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

    double MU = m_He * m_Ar / (m_He + m_Ar); 
    MoleculeSystem ms = init_ms(MU, ATOM, ATOM, NULL, NULL, seeds[_wrank]);

    dipole = dipole_lab;

    double tolerance = 1e-12;
    Trajectory traj = init_trajectory(&ms, tolerance);
   
    CalcParams params = {};
    params.ps                               = FREE_AND_METASTABLE;
    params.sampler_Rmin                     = 4.0;
    params.sampler_Rmax                     = 40.0;
    params.niterations                      = 100;
    params.total_trajectories               = 2000000;
    params.cvode_tolerance                  = 1e-12;
    params.sampling_time                    = 200.0;
    params.MaxTrajectoryLength              = 65536;
    params.Rcut                             = 40.0;
    params.partial_partition_function_ratio = 2.67619e+05;
    params.initialM0_npoints                = 10000000;
    params.initialM2_npoints                = 10000000;
    params.pesmin                            = V_HeAr(Rmin) / HTOCM;
    params.cf_filename                      = "./CF-He-Ar-F-300.0.txt";

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
