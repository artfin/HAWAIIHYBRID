#ifndef HAWAII_H_
#define HAWAII_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

#ifdef __cplusplus
extern "C" {
#endif

// This enum allows us to both differentiate between the systems of different type
// as well as to store the size of phase point: 
//   size(phase_point) = MonomerType % modulo_base
#define modulo_base 100
typedef enum {
    ATOM                                 = 0,
    LINEAR_MOLECULE                      = 4,
    LINEAR_MOLECULE_REQUANTIZED_ROTATION = 104,
    ROTOR                                = 6,
} MonomerType;

typedef struct { 
    MonomerType t;
    double I[3];
    double *qp;
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
    
    // TODO: I feel like it would be better to divide *dVdq array between intermolecular coordinates and monomers
    // just like qp
    double *dVdq;
    N_Vector dVdq_qp;
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

    double sampling_time;
    double R0; // initial distance for pr/mu calculation 
    double Rcut;
    size_t MaxTrajectoryLength;
    size_t CF_Length;
    // для расчета массива корреляционных функций
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


/* ----------------------------- */
/*    User-Supplied Functions    */
/* ------------------------------*/
typedef void (*dipolePtr)(double *q, double dip[3]);
extern dipolePtr dipole;

extern double pes(double *q);
extern void dpes(double *q, double *dVdq);
/* ------------------------------*/

MoleculeSystem init_ms(double mu, MonomerType t1, MonomerType t2, double *I1, double *I2, size_t seed); 
void free_ms(MoleculeSystem *ms);

const char* monomer_type_name(MonomerType t);

void put_qp_into_ms(MoleculeSystem *ms, Array qp);
void get_qp_from_ms(MoleculeSystem *ms, Array *qp);

// this function takes the phase point from MoleculeSystem 
// and packs it into "intermediate_q" field. Then this pointer
// can be passed to potential energy / its derivatives / dipole function 
void extract_q_and_write_into_ms(MoleculeSystem *ms);


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
double torque_monomer(MoleculeSystem *ms, size_t monomer_index); 

void calculate_M0(MoleculeSystem *ms, CalcParams *params, double Temperature, double *m, double *q);

int assert_float_is_equal_to(double estimate, double true_value, double abs_tolerance);



#ifdef __cplusplus
}
#endif

#endif // HAWAII_H_
