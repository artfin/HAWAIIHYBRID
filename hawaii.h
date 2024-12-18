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

#include "constants.h"

#define IPHI 0   
#define IPPHI 1  
#define ITHETA 2
#define IPTHETA 3
#define IR 4 
#define IPR 5    
// used only in Monomer.qp
#define IPSI 4
#define IPPSI 5

#define UNUSED(x) (void)(x)
#define TODO(message) do { fprintf(stderr, "%s:%d: TODO: %s\n", __FILE__, __LINE__, message); abort(); } while(0)
#define UNREACHABLE(message) do { fprintf(stderr, "%s:%d: UNREACHABLE: %s\n", __FILE__, __LINE__, message); abort(); } while(0)

typedef enum {
    ATOM = 0, 
    LINEAR_MOLECULE = 4, 
    ROTOR = 6,
} MonomerType;

typedef struct { 
    MonomerType t;
    double I[3];
    double *qp;
} Monomer; 

typedef struct {
    Monomer m1;
    Monomer m2;
    double mu;

    size_t Q_SIZE;
    size_t QP_SIZE;

    double intermolecular_qp[6];

    double *intermediate_q;
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

#include "array.h"

#ifdef __cplusplus
extern "C" {
#endif

extern double pes(double *q);
extern void dpes(double *q, double *dVdq);

void fill_qp(MoleculeSystem *ms, Array qp);

MoleculeSystem init_ms(double mu, MonomerType t1, MonomerType t2, double *I1, double *I2, size_t seed); 
void free_ms(MoleculeSystem *ms);

int rhs(realtype t, N_Vector y, N_Vector ydot, void *data);

double kinetic_energy(MoleculeSystem *ms);
double Hamiltonian(MoleculeSystem *ms);


double generate_normal(double sigma); 

void q_generator(MoleculeSystem *ms, CalcParams *params);
void p_generator(MoleculeSystem *ms, double Temperature);
bool reject(MoleculeSystem *ms, double Temperature, double pesmin);

typedef void (*dipolePtr)(double *q, double dip[3]);
extern dipolePtr dipole;

//extern void dipole(double *q, double dip[3]);
void calculate_M0(MoleculeSystem *ms, CalcParams *params, double Temperature, double *m, double *q);

#ifdef __cplusplus
}
#endif



#endif // HAWAII_H_
