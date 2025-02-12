#ifndef TRAJECTORY_H_
#define TRAJECTORY_H_

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "array.h"
#include "hawaii.h"

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <cvode/cvode_direct.h>
#include <sundials/sundials_types.h>


#ifdef __cplusplus
extern "C" {
#endif

extern int rhs(realtype t, N_Vector y, N_Vector ydot, void * user_data);

typedef struct
{
    size_t DIM;
    size_t mxsteps;
    N_Vector y;
    N_Vector abstol;
    realtype reltol; 
    SUNMatrix A;
    SUNLinearSolver LS;
    void *cvode_mem;
    //  
    MoleculeSystem* ms;
    N_Vector ic; 
    // -- energy conservation 
    bool check_energy;
    double E0;
    double E_last;
} Trajectory;

Trajectory init_trajectory(MoleculeSystem *ms, double reltol); 
void free_trajectory(Trajectory *traj);

int make_step(Trajectory *traj, double tout, double *t);

void set_initial_condition(Trajectory *traj, Array qp);
N_Vector make_vector(int size);
void set_tolerance(Trajectory *traj, double tolerance);

void reinit_trajectory(Trajectory *traj, double t);

#ifdef __cplusplus
}
#endif

#endif // TRAJECTORY_H_
