#include "trajectory.h"

Trajectory init_trajectory(size_t DIM, double reltol) 
{
    Trajectory t = {0};
    
    t.mxsteps = 50000;
    t.reltol = (realtype) reltol;

    t.DIM = DIM;
    t.y = make_vector(DIM);
    t.abstol = make_vector(DIM);

    // allocate solver memory and use Backward Differentiation formula
    t.cvode_mem = CVodeCreate(CV_BDF);
    assert(t.cvode_mem != NULL);

    int flag;

    // initialize the integrator memory and specify the right-hand side function, 
    // the initial time t0, and the initial dependent variable vector
    flag = CVodeInit(t.cvode_mem, rhs, 0.0 /* initial time */, t.y);
    if (flag < 0) {
        fprintf(stderr, "ERROR: CVodeInit -- unrecoverable error\n"); 
        exit(1);
    }

    t.A = SUNDenseMatrix(DIM, DIM);
    assert(t.A != NULL); 

    t.LS = SUNDenseLinearSolver(t.y, t.A);
    assert(t.LS != NULL);

    flag = CVDlsSetLinearSolver(t.cvode_mem, t.LS, t.A);
    if (flag < 0) {
        fprintf(stderr, "ERROR: CVDlsSetLinearSolver -- unrecoverable error\n"); 
        exit(1);
    }

    set_tolerance(&t, reltol);
    
    // number of internal steps before tout; default is 500
    CVodeSetMaxNumSteps(t.cvode_mem, t.mxsteps);

    // the function specifies the maximum number of messages issued by the solver warning that t + h = t on the next internal step
    // default value is 10. A negative value for mxhnil indicates that no warning messages should be issued 
    CVodeSetMaxHnilWarns(t.cvode_mem, -1);
}

void free_trajectory(Trajectory *t)
{
    N_VDestroy(t->y);
    N_VDestroy(t->abstol);

    CVodeFree(&t->cvode_mem);
    SUNLinSolFree(t->LS);
    SUNMatDestroy(t->A);
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

