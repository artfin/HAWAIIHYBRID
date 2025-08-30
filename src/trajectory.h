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
    size_t DIM; ///< Number of Hamilton's equations of motion (system dimension). 
    size_t mxsteps; ///< Maximum number of internal time steps for the solver.
    N_Vector y; ///< Internal phase-space vector (current step of the system).
    N_Vector abstol; ///< Absolute tolerance for each variable in the ODE system.
    realtype reltol; ///< Relative tolerance for the ODE solver.
    SUNMatrix A; ///< Matrix structure for storing jacobian information. 
    SUNLinearSolver LS; ///< Linear solver object for the ODE system.
    void *cvode_mem; ///< Pointer to the CVODE memory block (solver internal state).
    // 
    Array temp_qp;  ///< Temporary storage array for coordinates and momenta during trajectorypropagation.
    MoleculeSystem* ms; ///< Pointer to MoleculeSystem.
    N_Vector ic; ///< Initial conditions for the trajectory.
    // -- energy conservation 
    bool check_energy_conservation; ///< Flag indicating whether to monitor energy conservation during trajectory propagation.
    double E0; ///< Initial total energy of the system (used for conservation checks).
    double E_last; ///< Energy value from the last computed step (for tracking changes in energy).
} Trajectory;

Trajectory init_trajectory(MoleculeSystem *ms, double reltol); 
void free_trajectory(Trajectory *traj);

bool trajectory_apply_requantization(Trajectory *traj);
void trajectory_reinit(Trajectory *traj);
int make_step(Trajectory *traj, double tout, double *t);

void set_initial_condition(Trajectory *traj, Array qp);
N_Vector make_vector(int size);
void set_tolerance(Trajectory *traj, double tolerance);

void reinit_trajectory(Trajectory *traj, double t);

#ifdef __cplusplus
}
#endif

#endif // TRAJECTORY_H_
