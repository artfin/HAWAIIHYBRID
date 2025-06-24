#include "hawaii.h"

#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

#include "ai_pes_co2ar_lib.hpp"
#include "ai_ids_co2ar_lib.hpp"

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
        printf("dot(%s): \t %.10e \t %.10e \t %.10e\n", var_to_cstring(i), NV_Ith_S(ydot, i), num_derivatives.data[i], NV_Ith_S(ydot, i) - num_derivatives.data[i]);
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

void test_jac()
{
    printf("\n-----------------------------------------\n");
    printf("Testing analytic derivatives of angles transformation against the numerical ones\n");

    double qlab[5] = {0.1, 0.2, 5.0, 0.3, 0.4};

    size_t order = 6;
    size_t ninput_coordinates = 5;
    size_t noutput_coordinates = 5;
    gsl_matrix* numerical_jac = compute_numerical_jac(linear_molecule_atom_lab_to_mol, qlab, ninput_coordinates, noutput_coordinates, order);

    printf("Numerical jacobian:\n");
    for (size_t i = 0; i < ninput_coordinates; ++i) {
        for (size_t j = 0; j < noutput_coordinates; ++j) {
            printf("%.10e ", gsl_matrix_get(numerical_jac, i, j));
        }
        printf("\n");
    }
   
    static Eigen::Matrix<double, 5, 5> analytic_jac = Eigen::Matrix<double, 5, 5>::Zero(5, 5);
    double qmol[5];
    linear_molecule_atom_lab_to_mol(qlab, qmol);
    linear_molecule_atom_Jacobi_mol_by_lab(analytic_jac, qlab, qmol); 
    
    printf("Analytic jacobian:\n");
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            printf("%.10e ", analytic_jac(i, j)); 
        }
        printf("\n");
    }

    printf("Difference:\n");
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            printf("%.10e ", analytic_jac(i, j) - gsl_matrix_get(numerical_jac, i, j)); 
        }
        printf("\n");
    }

    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            if (assert_float_is_equal_to(analytic_jac(i, j), gsl_matrix_get(numerical_jac, i, j), 1e-9) > 0) {
                printf("ERROR: The elements (%zu, %zu) disagree!\n", i, j);
                exit(1);
            }
        }
    }

    printf("-----------------------------------------\n");

    gsl_matrix_free(numerical_jac);
}

int main()
{
    int seed = 42;

    double MU = m_CO2 * m_Ar / (m_CO2 + m_Ar); 
    double I1[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE, ATOM, I1, NULL, seed);

    pes_init();
    pes = pes_lab;
    dpes = dpes_lab;

    Array qp = create_array(ms.QP_SIZE);
    double data[] = {
        1.748408280024098e-01, // PHI
        8.000086892719947e+00, // PPHI
        2.043243374007463e-01, // THETA
        1.292311668112794e+01, // PTHETA
        19.066821669009366e+00, // R
       -0.186729452372973e-01, // PR
        1.462013919064840e+00, // PHI1
        0.199991310728005e+01, // PPHI1
        1.700437238862704e+00, // THETA1
        0.160851788770303e+01, // PTHETA1
    }; 
    
    init_array(&qp, data, ms.QP_SIZE);
    put_qp_into_ms(&ms, qp);

    test_jac();
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

    for (size_t nstep = 0; ; ++nstep, tout += params.sampling_time)
    {
        int status = make_step(&traj, tout, &t);

        if (status > 0) {
            fprintf(stderr, "CVODE ERROR: status = %d\n", status);
            exit(1);
        }
      
        double E = Hamiltonian(&ms);
        printf("%10.1lf \t %12.10lf \t %12.15lf\n", t, ms.intermolecular_qp[IR], E-E0);

        if (ms.intermolecular_qp[IR] > 40.0) break;
    }

    free_trajectory(&traj);
    free_ms(&ms);
    free_array(&qp);

    return 0;
}

