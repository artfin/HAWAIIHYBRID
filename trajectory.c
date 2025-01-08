#include "trajectory.h"

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

    int flag;
   
    Array qp = create_array(ms->QP_SIZE);
    get_qp_from_ms(ms, &qp);
    memcpy(N_VGetArrayPointer(traj.y), qp.data, ms->QP_SIZE * sizeof(double));
    free_array(&qp);

    // initialize the integrator memory and specify the right-hand side function, 
    // the initial time t0, and the initial dependent variable vector
    flag = CVodeInit(traj.cvode_mem, rhs, 0.0 /* initial time */, traj.y);
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

void free_trajectory(Trajectory *t)
{
    N_VDestroy(t->y);
    N_VDestroy(t->abstol);

    SUNLinSolFree(t->LS);
    SUNMatDestroy(t->A);

    CVodeFree(&t->cvode_mem);
}


int make_step(Trajectory *traj, double tout, double *t)
{
    int flag = CVode(traj->cvode_mem, tout, traj->y, t, CV_NORMAL);
    return flag;
}

void set_initial_condition(Trajectory *t, Array qp)
{
    memcpy(N_VGetArrayPointer(t->y), qp.data, qp.n * sizeof(double));
    CVodeReInit(t->cvode_mem, 0.0, t->y);
}

N_Vector make_vector(const int size) 
{
    N_Vector nv = N_VNew_Serial(size);
    assert(nv != NULL); 

    return nv;
}


void set_tolerance(Trajectory *t, double tolerance) 
{
    t->reltol = tolerance;
    for (size_t k = 0; k < t->DIM; ++k) { 
        NV_Ith_S(t->abstol, k) = tolerance; 
    }

    // specify the scalar relative tolerance and vector absolute tolerances
    int flag = CVodeSVtolerances(t->cvode_mem, t->reltol, t->abstol);
    if (flag < 0) {
        fprintf(stderr, "ERROR: CVodeSVtolerances -- unrecoverable error\n"); 
        exit(1);
    }
}

void reinit_trajectory(Trajectory *traj, double t) {
    CVodeReInit(traj->cvode_mem, t, traj->y);
}

