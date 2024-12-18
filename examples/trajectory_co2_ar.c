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


int main()
{
    int seed = 42;

    double MU = m_CO2 * m_Ar / (m_CO2 + m_Ar); 
    double I1[2] = {II_CO2, II_CO2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE, ATOM, I1, NULL, seed);

    co2_ar_pes.init();

    printf("QP_SIZE = %zu\n", ms.QP_SIZE);
    Array qp = create_array(ms.QP_SIZE);
    double data[] = {7.0, 8.0, 9.0, 10.0, 5.0, 6.0, 11.0, 12.0, 13.0, 14.0};
    init_array(&qp, data, ms.QP_SIZE);
    print_array(qp);

    fill_qp(&ms, qp);

    double reltol = 1e-15;
    Trajectory traj = init_trajectory(&ms, reltol);
   
    set_initial_condition(&traj, qp);

    CalcParams params = {};
    params.sampling_time = 200.0; 
    
    double t = 0.0;
    double tout = params.sampling_time;

    size_t nsteps = 100;

    for (size_t nstep = 0; nstep < nsteps; ++nstep, tout += params.sampling_time)
    {
        int status = make_step(&traj, tout, &t);
        if (status > 0) {
            fprintf(stderr, "CVODE ERROR: status = %d\n", status);
            exit(1);
        }
        
        init_array(&qp, N_VGetArrayPointer(traj.y), ms.QP_SIZE);
        fill_qp(&ms, qp);
        double E = Hamiltonian(&ms);
    
        printf("%.1lf %.10lf %.10lf\n", t, NV_Ith_S(traj.y, IR), E); 
    }

    // free_trajectory(&traj);

    return 0;
}

