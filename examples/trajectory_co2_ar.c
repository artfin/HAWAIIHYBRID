#include "hawaii.h"

#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

#include "ai_pes_co2ar.h"

double pes(double *q) {
    static double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);
    return pes_co2ar(qmol[0], qmol[4]);
} 

void dpes(double *q, double *dpesdq) {
    static Eigen::Matrix<double, 5, 5> jac;
    static Eigen::Matrix<double, 5, 1> derivatives_mol, derivatives_lab; 
    static double qmol[5];
    
    jac.setZero();
    linear_molecule_atom_lab_to_mol(q, qmol);
    linear_molecule_atom_Jacobi_mol_by_lab(jac, q, qmol);  
    
    double dR, dTheta;
    dpes_co2ar(qmol[0], qmol[4], &dR, &dTheta); // [R, THETAM] -> [dpes_dR, dpes_dTheta]

    // [dpes_dR, 0, 0, 0, dpes_dTheta]
    derivatives_mol(0) = dR; 
    derivatives_mol(4) = dTheta; 

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

    double MU = m_CO2 * m_Ar / (m_CO2 + m_Ar); 
    double I1[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE, ATOM, I1, NULL, seed);

    init_pes();

    Array qp = create_array(ms.QP_SIZE);
    double data[] = {0.1, 8.0, 0.2, 10.0, 5.0, 6.0, 0.3, 12.0, 0.4, 14.0};
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

    size_t nsteps = 200;

    for (size_t nstep = 0; nstep < nsteps; ++nstep, tout += params.sampling_time)
    {
        int status = make_step(&traj, tout, &t);
        if (status > 0) {
            fprintf(stderr, "CVODE ERROR: status = %d\n", status);
            exit(1);
        }
        
        double E = Hamiltonian(&ms);
        printf("%.1lf \t %.10lf \t %.15lf\n", t, NV_Ith_S(traj.y, IR), E-E0);

        if (assert_float_is_equal_to(E-E0, 0.0, 1e-8) > 0) {
            exit(1); 
        }
    }

    free_trajectory(&traj);
    free_ms(&ms);
    free_array(&qp);

    return 0;
}

