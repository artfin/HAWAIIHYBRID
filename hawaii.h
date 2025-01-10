#ifndef HAWAII_H_
#define HAWAII_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include <gsl/gsl_histogram.h>

#ifndef __cplusplus
#define _GNU_SOURCE
void sincos(double, double*, double*);
#endif
#include <math.h>

#ifndef __cplusplus
#define MT_GENERATE_CODE_IN_HEADER 0
#endif 
#include "mtwist.h"
/*
 * double mt_drand(void)
 *   Return a pseudorandom double in [0,1) with 32 bits of randomness
 *
 * uint32_t mt_lrand(void);
 *   Generate 32-bit random value 
 *
 */

#include <assert.h>
#include <stdbool.h>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>

#include "array.h"
#include "constants.h"

#define IPHI    0
#define IPPHI   1
#define ITHETA  2
#define IPTHETA 3
#define IR      4
#define IPR     5
// used only in Monomer.qp
#define IPSI    4
#define IPPSI   5

#define UNUSED(x) (void)(x)
#define TODO(message) do { fprintf(stderr, "%s:%d: TODO: %s\n", __FILE__, __LINE__, message); abort(); } while(0)
#define UNREACHABLE(message) do { fprintf(stderr, "%s:%d: UNREACHABLE: %s\n", __FILE__, __LINE__, message); abort(); } while(0)

#ifdef USE_MPI
#define INIT_WRANK                          \
    int _wrank;                             \
    MPI_Comm_rank(MPI_COMM_WORLD, &_wrank); \

#define INIT_WSIZE                          \
    int _wsize;                             \
    MPI_Comm_size(MPI_COMM_WORLD, &_wsize); \

#define PRINT0(...) if (_wrank == 0) printf(__VA_ARGS__); 
#else
#define INIT_WRANK
#define INIT_WSIZE
#define PRINT0(...) printf(__VA_ARGS__)
#endif

#ifdef __cplusplus
extern "C" {
#endif

// This enum allows us to both differentiate between the systems of different type
// as well as to store the size of phase point: 
//   size(phase_point) = MonomerType%MODULO_BASE
#define MODULO_BASE 100
typedef enum {
    ATOM                                 = 0,
    LINEAR_MOLECULE                      = 4,
    LINEAR_MOLECULE_REQUANTIZED_ROTATION = MODULO_BASE + 4,
    ROTOR                                = 6,
    ROTOR_REQUANTIZED_ROTATION           = MODULO_BASE + 6,
} MonomerType;

typedef struct { 
    MonomerType t;
    double I[3];
    double *qp;
    double *dVdq;
    // TODO: should this be a monomer field or a bool[2] in MoleculeSystem?
    bool apply_requantization;
} Monomer; 

typedef struct {
    double intermolecular_qp[6];
    Monomer m1;
    Monomer m2;
    double mu;

    size_t Q_SIZE;
    size_t QP_SIZE;

    double *intermediate_q;
    double *dVdq;
} MoleculeSystem;

typedef enum {
    FREE_AND_METASTABLE, // автоматически отделяем вклады
    BOUND,
} PairState;

typedef struct {
    PairState ps;

    /* sampling */
    double sampler_Rmin;
    double sampler_Rmax;
    double pesmin;

    /* initial spectral moments check */ 
    size_t initialM0_npoints;
    double partial_partition_function_ratio;
   
    /* requantization */ 
    size_t torque_cache_len;
    double torque_bound;

    /* trajectory */
    double sampling_time;
    size_t MaxTrajectoryLength;
    double cvode_tolerance;
    
    /* correlation & pr/mu calculation */
    const char* cf_filename;
    size_t total_trajectories;
    size_t niterations;
    double Rcut;
    double R0; // initial distance for pr/mu calculation 

    /* correlation function array */
    double *temperatures;
    size_t ntemperatures;
} CalcParams;


typedef enum {
    PRMU,
    CORRELATION_SINGLE,
    CORRELATION_ARRAY,
} CalculationType;

typedef enum {
    REJECT,
    MCMC
} SamplingType;


typedef struct {
    double *t; 
    double *data;
    size_t len;
    size_t ntraj; 
    double T;
} CFnc;

typedef struct {
    double before2;
    double before;
    double current;

    size_t turning_points;
    
    size_t called;
    bool empty;
} Tracker;


/* ----------------------------- */
/*    User-Supplied Functions    */
/* ------------------------------*/
typedef void (*dipolePtr)(double*, double[3]);
extern dipolePtr dipole;

extern double pes(double *q);
extern void dpes(double *q, double *dVdq);
/* ------------------------------*/

MoleculeSystem init_ms(double mu, MonomerType t1, MonomerType t2, double *I1, double *I2, size_t seed); 
void free_ms(MoleculeSystem *ms);

const char* monomer_type_name(MonomerType t);
const char* pair_state_name(PairState ps);

void put_qp_into_ms(MoleculeSystem *ms, Array qp);
void get_qp_from_ms(MoleculeSystem *ms, Array *qp);

// this function takes the phase point from MoleculeSystem 
// and packs it into "intermediate_q" field. Then this pointer
// can be passed to potential energy / its derivatives / dipole function 
void extract_q_and_write_into_ms(MoleculeSystem *ms);

// this function takes the corresponding components of "ms.dVdq" and writes
// them into "m.dVdq" vectors for each monomer 
void extract_dVdq_and_write_into_monomers(MoleculeSystem *ms);

// 20.12.2024 NOTE: 
// the MoleculeSystem struct needs to be passed as void* into "rhs" function.
// Then, the first step is to write the phase point (N_Vector y) into the fields of MoleculeSystem
// using the method "put_qp_into_ms".
int rhs(realtype t, N_Vector y, N_Vector ydot, void *data);
Array compute_numerical_rhs(MoleculeSystem *ms); 

double kinetic_energy(MoleculeSystem *ms);
double Hamiltonian(MoleculeSystem *ms);

double generate_normal(double sigma); 

void q_generator(MoleculeSystem *ms, CalcParams *params);
void p_generator(MoleculeSystem *ms, double Temperature);
bool reject(MoleculeSystem *ms, double Temperature, double pesmin);

double j_monomer(Monomer m);
double torque_monomer(Monomer m); 

void invert_momenta(MoleculeSystem *ms);

void calculate_M0(MoleculeSystem *ms, CalcParams *params, double Temperature, double *m, double *q);
double analytic_full_partition_function_by_V(MoleculeSystem *ms, double T);

#ifdef USE_MPI
void mpi_calculate_M0(MoleculeSystem *ms, CalcParams *params, double Temperature, double *m, double *q);
CFnc calculate_correlation_and_save(MoleculeSystem *ms, CalcParams *params, double Temperature);
#endif // USE_MPI
    
double* linspace(double start, double end, size_t n);
// std::vector<double> arange(double start, double step, size_t size);

void free_cfnc(CFnc cf);

void save_correlation_function(FILE *fd, CFnc crln, CalcParams *params);

int assert_float_is_equal_to(double estimate, double true_value, double abs_tolerance);

#ifdef __cplusplus
}
#endif

#endif // HAWAII_H_
