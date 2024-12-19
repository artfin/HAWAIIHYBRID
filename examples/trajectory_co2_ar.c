#include "hawaii.h"

#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

#include "ai_pes_co2ar.hpp"
static AI_PES_co2_ar co2_ar_pes;

double pes(double *q) {
    return 0.0;

    double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);
    return co2_ar_pes.pes(qmol[0], qmol[4]);
} 

void dpes(double *q, double *dq) {
    memset(dq, 0.0, 10 * sizeof(double)); 
    UNUSED(q);
    // TODO("dpes");
}

void test_rhs(MoleculeSystem *ms, Array qp)
{
    N_Vector y    = make_vector(ms->QP_SIZE);
    N_Vector ydot = make_vector(ms->QP_SIZE);

    memcpy(N_VGetArrayPointer(y), qp.data, qp.n * sizeof(double));
    rhs(0.0, y, ydot, (void*) ms);    

    Array num_derivatives = compute_numerical_rhs(ms);

    printf("# \t analytic \t numeric \t difference \n");
    for (size_t i = 0; i < ms->QP_SIZE; ++i) {
        printf("%zu: %.10e \t %.10e \t %.10e\n", i, NV_Ith_S(ydot, i), num_derivatives.data[i], NV_Ith_S(ydot, i) - num_derivatives.data[i]);

        if (assert_float_is_equal_to(NV_Ith_S(ydot, i), num_derivatives.data[i], 1e-9) > 0) {
            exit(1);
        }
    }
    printf("-----------------------------------------");
   
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

    co2_ar_pes.init();

    Array qp = create_array(ms.QP_SIZE);
    double data[] = {7.0, 8.0, 9.0, 10.0, -2.0, 6.0, 11.0, 12.0, 13.0, 14.0};
    init_array(&qp, data, ms.QP_SIZE);
    put_qp_into_ms(&ms, qp);

    test_rhs(&ms, qp);

    double reltol = 1e-12;
    Trajectory traj = init_trajectory(&ms, reltol);
   
    set_initial_condition(&traj, qp);

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
        
        init_array(&qp, N_VGetArrayPointer(traj.y), ms.QP_SIZE);
        put_qp_into_ms(&ms, qp);
        double E = Hamiltonian(&ms);
    
        printf("%.1lf %.10lf %.10lf\n", t, NV_Ith_S(traj.y, IR), E); 
    }

    free_trajectory(&traj);
    free_ms(&ms);
    free_array(&qp);

    return 0;
}

