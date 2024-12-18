#ifndef TRAJECTORY_H_
#define TRAJECTORY_H_

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "array.h"

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <cvode/cvode_direct.h>
#include <sundials/sundials_types.h>

int rhs(realtype t, N_Vector y, N_Vector ydot, void * user_data);

typedef struct
{
    size_t DIM;
    size_t mxsteps;
    N_Vector y         ;
    N_Vector abstol    ;
    realtype reltol    ; 
    SUNMatrix A        ;
    SUNLinearSolver LS ;
    void * cvode_mem   ;
} Trajectory;

#ifdef __cplusplus
extern "C" {
#endif

Trajectory init_trajectory(size_t DIM, double reltol);
void free_trajectory(Trajectory *t);

int make_step(Trajectory *traj, double tout, double *t);

void set_initial_condition(Trajectory *t, Array qp);
N_Vector make_vector(const int size);
void set_tolerance(Trajectory *t, double tolerance);

#ifdef __cplusplus
}
#endif

#endif // TRAJECTORY_H_
