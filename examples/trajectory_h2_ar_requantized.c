#include "hawaii.h"
#include "array.h"
#include "trajectory.h"
#include "angles_handler.hpp"

#include "ai_pes_h2ar_leroy.h"

#include <alloca.h>

double pes(double *q) {
    static double qmol[5];
    linear_molecule_atom_lab_to_mol(q, qmol);
    return pes_h2ar(qmol[0], qmol[4]);
} 

void dpes(double *q, double *dpesdq) {
    static Eigen::Matrix<double, 5, 5> jac;
    static Eigen::Matrix<double, 5, 1> derivatives_mol, derivatives_lab; 
    static double qmol[5];
    
    jac.setZero();
    linear_molecule_atom_lab_to_mol(q, qmol);
    linear_molecule_atom_Jacobi_mol_by_lab(jac, q, qmol);  
    
    double dR, dTheta;
    dpes_h2ar(qmol[0], qmol[4], &dR, &dTheta); // [R, THETAM] -> [dpes_dR, dpes_dTheta]

    // [dpes_dR, 0, 0, 0, dpes_dTheta]
    derivatives_mol(0) = dR; 
    derivatives_mol(4) = dTheta; 
    
    derivatives_lab = jac * derivatives_mol;
    Eigen::VectorXd::Map(dpesdq, 5) = derivatives_lab;
}

#define IPHI1 6
#define IPPHI1 7
#define ITHETA1 8
#define IPTHETA1 9

int main()
{
    int seed = 42;

    double MU = m_H2 * m_Ar / (m_H2 + m_Ar); 
    double I1[2] = {II_H2, II_H2};
    MoleculeSystem ms = init_ms(MU, LINEAR_MOLECULE_REQUANTIZED_ROTATION, ATOM, I1, NULL, seed);

    Array qp = create_array(ms.QP_SIZE);
    double data[ms.QP_SIZE];
    data[IR]       = 10.0000000000000000;
    data[IPR]      = -0.1390728108333908;
    data[IPHI]     = 2.6149492227311319;
    data[IPPHI]    = -0.845274877844076;
    data[ITHETA]   = 1.5712832054878647;
    data[IPTHETA]  = 0.7051209557879594;
    data[IPHI1]    = 1.7133086475792336;
    data[IPPHI1]   = -0.4428409546240073;
    data[ITHETA1]  = 0.5404197109595652;
    data[IPTHETA1] = 0.3298923262808994;
    
    init_array(&qp, data, ms.QP_SIZE);
    put_qp_into_ms(&ms, qp);

    double tolerance = 1e-15;
    Trajectory traj = init_trajectory(&ms, tolerance);
   
    double E0 = Hamiltonian(&ms);
    set_initial_condition(&traj, qp);
    print_array(qp);

    CalcParams params = {};
    params.sampling_time = 200.0; 
    params.torque_cache_len = 10;
    params.torque_bound = 1e-5;
    
    size_t switch_counter = 0;

    double t = 0.0;
    double tout = params.sampling_time;

    double* cache = (double*) alloca(params.torque_cache_len * sizeof(double));
    size_t cur_cache = 0;

    for (size_t nstep = 0; ; ++nstep, tout += params.sampling_time)
    {
        int status = make_step(&traj, tout, &t);
        if (status > 0) {
            fprintf(stderr, "CVODE ERROR: status = %d\n", status);
            exit(1);
        }
        
        put_qp_into_ms(&ms, (Array){.data = N_VGetArrayPointer(traj.y), .n = ms.QP_SIZE});
    
        double E = Hamiltonian(&ms);
        double j = j_monomer(ms.m1);
        double torq = torque_monomer(ms.m1); 

        cache[nstep % params.torque_cache_len] = torq;
        bool all_less_than_bound = true;
        bool all_more_than_bound = true;
        for (size_t i = 0; i < params.torque_cache_len; ++i) {
            if (torq > params.torque_bound) {
                all_less_than_bound = false;
            }
            if (torq < params.torque_bound) {
                all_more_than_bound = false;
            } 
        }

        if (all_less_than_bound) {
            if (!ms.m1.apply_requantization) {
                ms.m1.apply_requantization = true;
                switch_counter++;

                printf("Setting requantization to 'true': switch counter = %zu\n", switch_counter);
            }
        }

        if (all_more_than_bound) {
            if (ms.m1.apply_requantization) {
                ms.m1.apply_requantization = false;
                switch_counter++;
                
                printf("Setting requantization to 'false': switch counter = %zu\n", switch_counter);
            }
        }

        if (ms.m1.apply_requantization) {
            for (size_t i = 0; i < 4; ++i) {
                NV_Ith_S(traj.y, i + 6) = ms.m1.qp[i];  
            } 

            reinit_trajectory(&traj, tout);
        }


        printf("%10.1lf \t %12.10lf \t %12.15lf \t %12.5e \t %12.5e\n", t, ms.intermolecular_qp[IR], E-E0, j, torq);

        if (ms.intermolecular_qp[IR] > 40.0) break;
    }

    free_trajectory(&traj);
    free_ms(&ms);
    free_array(&qp);
    
    return 0;
}

