#include "trajectory.h"

#define TOLERANCE_FACTOR 1000 
   
Trajectory init_trajectory(MoleculeSystem *ms, double reltol) 
{
    Trajectory traj = {0};
    
    traj.mxsteps = 50000;
    traj.reltol = (realtype) reltol;

    traj.DIM = ms->QP_SIZE;
    traj.y = make_vector(ms->QP_SIZE);
    traj.abstol = make_vector(ms->QP_SIZE);

    // allocate solver memory and use Backward Differentiation formula
    traj.cvode_mem = CVodeCreate(CV_BDF);
    assert(traj.cvode_mem != NULL);

    traj.ms = ms;

    traj.temp_qp = create_array(ms->QP_SIZE);
    get_qp_from_ms(ms, &traj.temp_qp);
    memcpy(N_VGetArrayPointer(traj.y), traj.temp_qp.data, ms->QP_SIZE * sizeof(double));

    traj.check_energy_conservation = true;

    traj.ic = make_vector(ms->QP_SIZE);
    for (size_t i = 0; i < ms->QP_SIZE; ++i) {
        NV_Ith_S(traj.ic, i) = traj.temp_qp.data[i];
    } 

    // initialize the integrator memory and specify the right-hand side function, 
    // the initial time t0, and the initial dependent variable vector
    int flag = CVodeInit(traj.cvode_mem, rhs, 0.0 /* initial time */, traj.y);
    if (flag < 0) {
        fprintf(stderr, "SUNDIALS ERROR: CVodeInit -- unrecoverable error\n"); 
        exit(1);
    }

    traj.A = SUNDenseMatrix(ms->QP_SIZE, ms->QP_SIZE);
    assert(traj.A != NULL); 

    traj.LS = SUNDenseLinearSolver(traj.y, traj.A);
    assert(traj.LS != NULL);

    flag = CVDlsSetLinearSolver(traj.cvode_mem, traj.LS, traj.A);
    if (flag < 0) {
        fprintf(stderr, "SUNDIALS ERROR: CVDlsSetLinearSolver -- unrecoverable error\n"); 
        exit(1);
    }

    set_tolerance(&traj, reltol);
    
    // number of internal steps before tout; default is 500
    CVodeSetMaxNumSteps(traj.cvode_mem, traj.mxsteps);

    // the function specifies the maximum number of messages issued by the solver warning that t + h = t on the next internal step
    // default value is 10. A negative value for mxhnil indicates that no warning messages should be issued 
    CVodeSetMaxHnilWarns(traj.cvode_mem, -1);

    flag = CVodeSetUserData(traj.cvode_mem, (void*) ms);
    if (flag < 0) {
        fprintf(stderr, "SUNDIALS ERROR: CVodeSetUserData -- unrecoverable error\n");
        exit(1);
    }

    return traj;
}

void free_trajectory(Trajectory *traj)
{
    N_VDestroy(traj->y);
    N_VDestroy(traj->abstol);
    N_VDestroy(traj->ic);

    SUNLinSolFree(traj->LS);
    SUNMatDestroy(traj->A);

    CVodeFree(&traj->cvode_mem);
    free_array(&traj->temp_qp);
}


int make_step(Trajectory *traj, double tout, double *t)
{
    int flag = CVode(traj->cvode_mem, tout, traj->y, t, CV_NORMAL);
    
    // 21.12.2024 NOTE: 
    // We copy the "N_Vector y" from "Trajectory" into "MoleculeSystem" on each call of rhs.
    // However after "rhs" has been executed, CVode advances the dynamic variables, therefore we 
    // need to update the phase-point in MoleculeSystem once CVode function completes  
    put_qp_into_ms(traj->ms, (Array){.data = N_VGetArrayPointer(traj->y), .n = traj->ms->QP_SIZE});

    if (traj->ms->m1.apply_requantization) {
        if (traj->ms->m1.t == LINEAR_MOLECULE_REQ_INTEGER) {
            double j[3];
            Monomer *m = &traj->ms->m1;

            j_monomer(*m, j);
            double jl   = sqrt(j[0]*j[0] + j[1]*j[1] + j[2]*j[2]);
            double jreq = find_closest_integer(jl);

            double scaling_factor = 0.0;
            if (jl > 1e-15) {
                scaling_factor = jreq / jl; 
            }

            m->qp[IPPHI]   *= scaling_factor;
            m->qp[IPTHETA] *= scaling_factor;

            //printf("make_step: pphi = %.5e => j after scaling = %.5e\n", m->qp[IPPHI], j_monomer(*m));
        }
        if (traj->ms->m1.t == LINEAR_MOLECULE_REQ_HALFINTEGER) {
            double j[3];
            Monomer *m = &traj->ms->m1;

            j_monomer(*m, j);
            double jl   = sqrt(j[0]*j[0] + j[1]*j[1] + j[2]*j[2]);
            double jreq = find_closest_half_integer(jl);

            double scaling_factor = 0.0;
            if (jl > 1e-15) {
                scaling_factor = jreq / jl; 
            }

            m->qp[IPPHI]   *= scaling_factor;
            m->qp[IPTHETA] *= scaling_factor;

            //printf("make_step: pphi = %.5e => j after scaling = %.5e\n", m->qp[IPPHI], j_monomer(*m));
        }

        get_qp_from_ms(traj->ms, &traj->temp_qp);
        memcpy(N_VGetArrayPointer(traj->y), traj->temp_qp.data, traj->ms->QP_SIZE * sizeof(double));

        CVodeReInit(traj->cvode_mem, *t, traj->y);
    }



    if (traj->check_energy_conservation) {
        traj->E_last = Hamiltonian(traj->ms);

        double abserr = fabs(traj->E_last - traj->E0);
        if (abserr > TOLERANCE_FACTOR * traj->reltol) {
            printf("WARNING: the change of energy (%.5e) during trajectory propagation exceeded desired tolerance (%.5e)\n",
                    abserr, traj->reltol);
            printf("The starting initial condition is: ");
            for (size_t i = 0; i < traj->DIM; ++i) {
                printf("%.10e ", NV_Ith_S(traj->ic, i));
            } 
            printf("\n");
        } 
    }

    return flag;
}

void set_initial_condition(Trajectory *traj, Array qp)
{
    memcpy(N_VGetArrayPointer(traj->y), qp.data, qp.n * sizeof(double));
    CVodeReInit(traj->cvode_mem, 0.0, traj->y);
    
    assert(qp.n == traj->DIM); 
    for (size_t i = 0; i < qp.n; ++i) {
        NV_Ith_S(traj->ic, i) = qp.data[i];
    } 
    
    put_qp_into_ms(traj->ms, (Array){.data = N_VGetArrayPointer(traj->y), .n = traj->ms->QP_SIZE});
    traj->E0 = Hamiltonian(traj->ms);
}

N_Vector make_vector(int size) 
{
    N_Vector nv = N_VNew_Serial(size);
    assert(nv != NULL); 

    return nv;
}


void set_tolerance(Trajectory *traj, double tolerance) 
{
    traj->reltol = tolerance;
    for (size_t k = 0; k < traj->DIM; ++k) { 
        NV_Ith_S(traj->abstol, k) = tolerance; 
    }

    // specify the scalar relative tolerance and vector absolute tolerances
    int flag = CVodeSVtolerances(traj->cvode_mem, traj->reltol, traj->abstol);
    if (flag < 0) {
        fprintf(stderr, "ERROR: CVodeSVtolerances -- unrecoverable error\n"); 
        exit(1);
    }
}

void reinit_trajectory(Trajectory *traj, double t) {
    CVodeReInit(traj->cvode_mem, t, traj->y);
}

