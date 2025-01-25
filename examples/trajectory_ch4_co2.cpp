#include "hawaii.h"

#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

#include "ai_pes_ch4_co2.h"
#include "ai_ids_ch4_co2.hpp"

#define NLAB 8
#define NKAL 5

double pes(double *qlab) {
    static double qkal[NKAL];
    CH4_linear_molecule_lab_to_kal(qlab, qkal); 
    return pes_ch4co2(qkal[0], qkal[1], qkal[2], qkal[3], qkal[4]); 
} 

void dpes(double *qlab, double *dpesdq) 
{
    static Eigen::Matrix<double, NLAB, NKAL> jac;
    static Eigen::Matrix<double, NKAL, 1> derivatives_kal;
    static Eigen::Matrix<double, NLAB, 1> derivatives_lab; 
    static double qkal[NKAL];

    jac.setZero();
    CH4_linear_molecule_lab_to_kal(qlab, qkal);
    CH4_linear_molecule_Jacobi_kal_by_lab(jac, qlab, qkal);  

    double dR, dphi1, dphi2, dtheta1, dtheta2; 
    dpes_ch4co2(qkal[0], qkal[1], qkal[2], qkal[3], qkal[4], &dR, &dphi1, &dtheta1, &dphi2, &dtheta2); 

    derivatives_kal(0) = dR; 
    derivatives_kal(1) = dphi1;
    derivatives_kal(2) = dtheta1;
    derivatives_kal(3) = dphi2;
    derivatives_kal(4) = dtheta2;

    derivatives_lab = jac * derivatives_kal;

    Eigen::VectorXd::Map(dpesdq, NLAB) = derivatives_lab;
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

const char *var_to_cstring(int n) {
    switch (n) {
        case 0:  return "Phi";
        case 1:  return "pPhi";
        case 2:  return "Theta";
        case 3:  return "pTheta";
        case 4:  return "R";
        case 5:  return "pR";
        case 6:  return "phi1t";
        case 7:  return "pPhi1t";
        case 8:  return "theta1t";
        case 9:  return "pTheta1t";
        case 10: return "psi1t";
        case 11: return "pPsi1t";
        case 12: return "phi2t";
        case 13: return "pPhi2t";
        case 14: return "theta2t";
        case 15: return "pTheta2t";
    }

    UNREACHABLE("var_to_cstring");
}

void test_jac()
{
    double qlab[8] = {0.1, 1.2, 8.0, 0.3, 0.4, 1.5, 0.6, 2.7};

    printf("\n-----------------------------------------\n");
    printf("Testing analytic derivatives of angles transformation against the numerical ones\n");

    size_t order = 6;
    size_t ninput_coordinates = 8;
    size_t noutput_coordinates = 5;
    gsl_matrix* numerical_jac = compute_numerical_jac(CH4_linear_molecule_lab_to_kal, qlab, ninput_coordinates, noutput_coordinates, order);

    printf("Numerical jacobian:\n");
    for (size_t i = 0; i < 8; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            printf("%.10e ", gsl_matrix_get(numerical_jac, i, j));
        }
        printf("\n");
    }
   
    static Eigen::Matrix<double, 8, 5> analytic_jac = Eigen::Matrix<double, 8, 5>::Zero(8, 5);
    double qkal[5];
    CH4_linear_molecule_lab_to_kal(qlab, qkal);
    CH4_linear_molecule_Jacobi_kal_by_lab(analytic_jac, qlab, qkal); 
    
    printf("Analytic jacobian:\n");
    for (size_t i = 0; i < 8; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            printf("%.10e ", analytic_jac(i, j)); 
        }
        printf("\n");
    }

    printf("Difference:\n");
    for (size_t i = 0; i < 8; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            printf("%.10e ", analytic_jac(i, j) - gsl_matrix_get(numerical_jac, i, j)); 
        }
        printf("\n");
    }

    for (size_t i = 0; i < 8; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            if (assert_float_is_equal_to(analytic_jac(i, j), gsl_matrix_get(numerical_jac, i, j), 1e-9) > 0) {
                printf("ERROR: The elements (%zu, %zu) disagree!\n", i, j);
                exit(1);
            }
        }
    }

    printf("-----------------------------------------\n");
    
    double dR, dphi1, dphi2, dtheta1, dtheta2; 
    dpes_ch4co2(qkal[0], qkal[1], qkal[2], qkal[3], qkal[4], &dR, &dphi1, &dtheta1, &dphi2, &dtheta2); 
    
    printf("\n-----------------------------------------\n");
    printf("Value of potential: %.10e\n", pes(qlab));
    printf("Value of derivatives in kalugina frame of reference:\n");
    printf("dV/dR:      %.10e\n", dR);
    printf("dV/dphi1:   %.10e\n", dphi1);
    printf("dV/dtheta1: %.10e\n", dtheta1);
    printf("dV/dphi2:   %.10e\n", dphi2);
    printf("dV/dtheta2: %.10e\n", dtheta2);
    printf("-----------------------------------------\n\n");
    
    Eigen::VectorXd derivatives_kal = Eigen::VectorXd::Zero(5, 1);
    derivatives_kal(0) = dR; 
    derivatives_kal(1) = dphi1;
    derivatives_kal(2) = dtheta1;
    derivatives_kal(3) = dphi2;
    derivatives_kal(4) = dtheta2;

    Eigen::VectorXd analytic_derivatives = analytic_jac * derivatives_kal;
   
    Array numeric_derivatives = compute_numerical_derivatives(pes, qlab, 8, 6);

    printf("-----------------------------------------\n");
    printf("Numeric vs analytic derivatives of potential energy in laboratory frame:\n");
    printf("Numeric \t analytic \t difference\n"); 
    for (size_t i = 0; i < 8; ++i) {
        double num = numeric_derivatives.data[i];
        double an = analytic_derivatives(i); 
        printf("%16.10e \t %16.10e \t %16.10e\n", num, an, num - an);
    } 
    printf("-----------------------------------------\n");

    free_array(&numeric_derivatives);
    gsl_matrix_free(numerical_jac);
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
        printf("dot(%s): \t %.10e \t %.10e \t %.10e\n", var_to_cstring(i), an, num, an - num);
    }
    
    for (size_t i = 0; i < ms->QP_SIZE; ++i) {
        if (assert_float_is_equal_to(NV_Ith_S(ydot, i), num_derivatives.data[i], 1e-4) > 0) {
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


int main()
{
    int seed = 42;

    double MU = m_CH4 * m_CO2 / (m_CH4 + m_CO2); 
    double I1[3] = {II_CH4, II_CH4, II_CH4};
    double I2[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(MU, ROTOR, LINEAR_MOLECULE, I1, I2, seed);

    init_pes();
    init_ids();
    
    Array qp = create_array(ms.QP_SIZE);
    
    double data[ms.QP_SIZE] = {
      0.1, -0.0115, 1.2, -0.037, 8.0, -0.03,   
      0.3, -2.622135, 0.4, -2.55059, 1.5, 4.0234,
      0.6, -0.033, 2.7, 0.51,
    }; 
    
    init_array(&qp, data, ms.QP_SIZE);
    put_qp_into_ms(&ms, qp);

    test_jac();
    test_rhs(&ms, qp);
        
    double tolerance = 1e-15;
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
        
        printf("%10.1lf \t %12.10lf \t %12.15lf %12.12lf\n", t, ms.intermolecular_qp[IR], E-E0, dipt[0]*dipt[0] + dipt[1]*dipt[1] + dipt[2]*dipt[2]);

        if (ms.intermolecular_qp[IR] > 30.0) break;

        //if (assert_float_is_equal_to(E-E0, 0.0, 1e-8) > 0) {
        //    exit(1); 
        //}
    }

    free_trajectory(&traj);
    free_ms(&ms);
    free_array(&qp);
    free_ids();

    return 0;
}

