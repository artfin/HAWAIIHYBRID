#include "hawaii.h"

#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

#include "ai_pes_ch4_co2.h"

double pes(double *qlab) {
    static double qkal[5];
    CH4_linear_molecule_lab_to_kal(qlab, qkal); 
    return pes_ch4co2(qkal[0], qkal[1], qkal[2], qkal[3], qkal[4]); 
} 

void dpes(double *q, double *dpesdq) {
    UNUSED(q);
    memset(dpesdq, 0, 8 * sizeof(double)); 
    // static Eigen::Matrix<double, 5, 5> jac;
    // static Eigen::Matrix<double, 5, 1> derivatives_mol, derivatives_lab; 
    // static double qmol[5];
    // 
    // jac.setZero();
    // linear_molecule_atom_lab_to_mol(q, qmol);
    // linear_molecule_atom_Jacobi_mol_by_lab(jac, q, qmol);  
    // 
    // double dR, dTheta;
    // dpes_co2ar(qmol[0], qmol[4], &dR, &dTheta); // [R, THETAM] -> [dpes_dR, dpes_dTheta]

    // // [dpes_dR, 0, 0, 0, dpes_dTheta]
    // derivatives_mol(0) = dR; 
    // derivatives_mol(4) = dTheta; 

    // derivatives_lab = jac * derivatives_mol;
    // Eigen::VectorXd::Map(dpesdq, 5) = derivatives_lab;
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

    Array num_derivatives = compute_numerical_rhs(ms);

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


int main()
{
    int seed = 42;

    double MU = m_CH4 * m_CO2 / (m_CH4 + m_CO2); 
    double I1[3] = {II_CH4, II_CH4, II_CH4};
    double I2[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(MU, ROTOR, LINEAR_MOLECULE, I1, I2, seed);

    init_pes();

    Array qp = create_array(ms.QP_SIZE);

    double data[ms.QP_SIZE] = {
      0.33, -0.0115, 0.44, -0.037, 25.0, -0.03,   
      0.55, -2.622135, 0.66, -2.55059, 0.77, 4.0234,
      0.88, -0.033, 0.99, 0.51,
    }; 
    
    init_array(&qp, data, ms.QP_SIZE);
    put_qp_into_ms(&ms, qp);
    
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
        printf("%10.1lf \t %12.10lf \t %12.15lf\n", t, ms.intermolecular_qp[IR], E-E0);

        if (ms.intermolecular_qp[IR] > 30.0) break;

        //if (assert_float_is_equal_to(E-E0, 0.0, 1e-8) > 0) {
        //    exit(1); 
        //}
    }

    free_trajectory(&traj);
    free_ms(&ms);
    free_array(&qp);

    return 0;
}

