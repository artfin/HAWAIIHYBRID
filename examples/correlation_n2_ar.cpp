#define USE_MPI
#include "hawaii.h"

#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

#include "ai_pes_n2_ar_pip_nn.hpp"
#include "ai_ids_n2_ar_pip_nn.hpp"

double pes(double *q)
// PHI THETA R PHI1T THETA1T  
{
    XYZ xyz = lab_to_xyz(q);
    return pes_xyz(&xyz); 
} 

void dpes(double *q, double *dpesdq) 
// PHI THETA R PHI1T THETA1T  
{
    XYZ xyz = lab_to_xyz(q);    
    XYZ dxyz{}; 

    dpes_xyz(&xyz, &dxyz);
    
    double Phi     = q[0];
    double Theta   = q[1];
    double R       = q[2];
    double phi1t   = q[3];
    double theta1t = q[4];
    
    double sinPhi   = sin(Phi);
    double cosPhi   = cos(Phi);
    double sinTheta = sin(Theta);
    double cosTheta = cos(Theta);
    
    /* dVdPhi */   dpesdq[0] = R * (-sinPhi * sinTheta * dxyz.c[6] + cosPhi * sinTheta * dxyz.c[7]) / HTOCM;
    /* dVdTheta */ dpesdq[1] = R * ( cosPhi * cosTheta * dxyz.c[6] + sinPhi * cosTheta * dxyz.c[7] - sinTheta * dxyz.c[8]) / HTOCM;
    /* dVdR */     dpesdq[2] = (cosPhi * sinTheta * dxyz.c[6] + sinPhi * sinTheta * dxyz.c[7] + cosTheta * dxyz.c[8]) / HTOCM;
    
    double sinPhi1t = sin(phi1t);
    double cosPhi1t = cos(phi1t);
    double sinTheta1t = sin(theta1t);
    double cosTheta1t = cos(theta1t);
    
    /* dVdphi1t */   dpesdq[3] = L_N2/2.0 * (-sinPhi1t * sinTheta1t * dxyz.c[0] + cosPhi1t * sinTheta1t * dxyz.c[1] + \
                                              sinPhi1t * sinTheta1t * dxyz.c[3] - cosPhi1t * sinTheta1t * dxyz.c[4]) / HTOCM;
    /* dVdtheta1t */ dpesdq[4] = L_N2/2.0 * (cosPhi1t * cosTheta1t * dxyz.c[0] + sinPhi1t * cosTheta1t * dxyz.c[1] - sinTheta1t * dxyz.c[2] - \
                                             cosPhi1t * cosTheta1t * dxyz.c[3] - sinPhi1t * cosTheta1t * dxyz.c[4] + sinTheta1t * dxyz.c[5]) / HTOCM;
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

void test_derivatives() 
{
    pes_init(true);

    double q[] = {0.1, 0.2, 8.0, 0.3, 0.4};
    double pv = pes(q);  
    printf("pes value = %.10e\n", pv);

    double dpesdq[5];
    dpes(q, dpesdq);

    for (size_t i = 0; i < 5; ++i) {
        double qq[5];
        memcpy(qq, q, 5 * sizeof(double));

        double delta = 1e-6;
        qq[i] = q[i] + delta; 
        double dp = pes(qq);

        qq[i] = q[i] - delta; 
        double dm = pes(qq);

        double num = (dp - dm) / (2.0*delta) / HTOCM; 
        printf("%zu: num = %.10e, analytical = %.10e\n", i, num, dpesdq[i]);
    }
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

    printf("#       \t analytic \t numeric \t difference \n");
    for (size_t i = 0; i < ms->QP_SIZE; ++i) {
        double an = NV_Ith_S(ydot, i);
        double num = num_derivatives.data[i]; 
        printf("%zu: \t %.10e \t %.10e \t %.10e\n", i, an, num, an - num);
    }
    
    for (size_t i = 0; i < ms->QP_SIZE; ++i) {
        if (assert_float_is_equal_to(NV_Ith_S(ydot, i), num_derivatives.data[i], 1e-3) > 0) {
            printf("ERROR: The element (%zu) disagree!\n", i);
            exit(1);
        }
    }
    printf("-----------------------------------------\n\n\n");
    
    put_qp_into_ms(ms, qp);
    
    free_array(&num_derivatives); 
    N_VDestroy(y);
    N_VDestroy(ydot);
}

void run_trajectory()
{
    double MU = m_N2 * m_Ar / (m_N2 + m_Ar); 
    double I1[2] = {II_N2, II_N2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE, ATOM, I1, NULL, 0);
    
    pes_init(true); 
    dipole_init(true);

    Array qp = create_array(ms.QP_SIZE);
    
    double data[ms.QP_SIZE] = {
      0.1, -0.0115, 1.2, -0.037, 13.0, -0.03,   
      0.3, -2.622135, 0.4, -2.55059, 
    }; 
    
    init_array(&qp, data, ms.QP_SIZE);
    put_qp_into_ms(&ms, qp);

    test_rhs(&ms, qp);
   
    double tolerance = 1e-12;
    Trajectory traj = init_trajectory(&ms, tolerance);
    
    double E0 = Hamiltonian(&ms);
    set_initial_condition(&traj, qp);
    
    double sampling_time = 200.0; 
    
    double t = 0.0;
    double tout = sampling_time;

    // 01.05 dipole values are the same as resulting from Supplemetary material of 
    // Continuum absorption in pure N2 gas and in its mixture with Ar, 10.1016/j.jqsrt.2024.109172    
    double dipt[3];
   
    for (size_t nstep = 0; ; ++nstep, tout += sampling_time)
    {
        int status = make_step(&traj, tout, &t);

        if (status > 0) {
            fprintf(stderr, "CVODE ERROR: status = %d\n", status);
            exit(1);
        }
        
        double E = Hamiltonian(&ms);

        extract_q_and_write_into_ms(&ms);
        dipole_lab(ms.intermediate_q, dipt);
        
        printf("%10.1lf \t %12.10lf \t %12.15lf \t %12.15lf\n", t, ms.intermolecular_qp[IR], E-E0, 
                dipt[0]*dipt[0]+dipt[1]*dipt[1]+dipt[2]*dipt[2]);

        if (ms.intermolecular_qp[IR] > 30.0) break;
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    INIT_WRANK;
    INIT_WSIZE;
    
    int seeds[] = {20106943, 74537365, 68838871, 45228894, 87924815, 75409565, 78079463, 71551420, 61270358, 12907113, 37812429, 43846516, 78563783, 67020198, 46463763, 15913615, 56260494, 73487105, 39711609, 17700181, 26290757, 10277061, 79863361, 29916933, 49969708, 56772087, 50637984, 10326549, 38641292, 81352997, 39721216, 73377986, 25410770, 65662596, 95180457, 12378249, 90265129, 87941783, 28243843, 40892493, 60538578, 67701455, 54511007, 95983369, 82062787, 59498255, 74798224, 29737064, 62225324, 18939487, 99787768, 87346783, 85009677, 63086756, 14433074, 99516596, 75192796, 75650862, 35503295, 20831228, 33162524, 12437412, 91210209, 76571171, 17922249, 15052896, 20751049, 49130285, 47619032, 29471832, 37775141, 74794610, 29101489, 55424490, 18716635, 20356657, 32172602, 97293267, 77861181, 29092863, 62639840, 56987692, 81112748, 58656852, 27469169, 27154689, 49220900, 95017234, 62779012, 97324386, 96580635, 45610890, 41955768, 14529933, 46207757, 44406420, 65851255, 34626188, 40063013, 66832671};
    assert(_wsize <= (int) (sizeof(seeds)/sizeof(int)) && "ASSERT: more seeds should be added\n");

    double MU = m_N2 * m_Ar / (m_N2 + m_Ar); 
    double I1[2] = {II_N2, II_N2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE, ATOM, I1, NULL, seeds[_wrank]);

    pes_init(true); 
    dipole_init(true);
    dipole = dipole_lab;

    double tolerance = 1e-12;
    Trajectory traj = init_trajectory(&ms, tolerance);

    CalcParams params = {};
    params.ps                               = FREE_AND_METASTABLE;
    params.sampler_Rmin                     = 4.51;
    params.sampler_Rmax                     = 40.0;
    params.niterations                      = 100;
    params.total_trajectories               = 10000000;
    params.cvode_tolerance                  = 1e-12;
    params.sampling_time                    = 200.0;
    params.MaxTrajectoryLength              = 65536;
    params.Rcut                             = 40.0;
    params.partial_partition_function_ratio = 2.68234e+05;
    params.initialM0_npoints                = 30000000;
    params.initialM2_npoints                = 30000000;
    params.pesmin                           = -97.3 / HTOCM;
    params.cf_filename                      = "./CF-N2-Ar-F-300.0.txt";

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


}
