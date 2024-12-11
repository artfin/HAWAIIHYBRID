#ifndef HAWAII_H_
#define HAWAII_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>

#define IPHI 0   
#define IPPHI 1  
#define ITHETA 2
#define IPTHETA 3
#define IR 4 
#define IPR 5    

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

typedef struct {
    double sampling_time;
    double R0;
    double Rcut;
    size_t MaxTrajectoryLength;
    size_t CF_Length;
    // для расчета массива корреляционных функций
    double *temperatures;
    size_t ntemperatures;
} CalcParams;

typedef enum {
    FREE_AND_METASTABLE, // автоматически отделяем вклады
    BOUND,
} PairState;

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
    double *data;
    size_t n;
} Array;

#define RAMTOAMU 1822.888485332
    
#define m_C  12.000000000000 * RAMTOAMU
#define m_N  14.003074004460 * RAMTOAMU
#define m_O  15.994914619598 * RAMTOAMU
#define m_Ar 39.9623831237 * RAMTOAMU

#define m_CO2 (m_C + 2.0 * m_O)
#define l_CO2 4.398

#define II_CO2  m_O / 2.0 * l_CO2 * l_CO2


extern double pes(double *q);
extern void dpes(double *q, double *dVdq);

Array create_array(size_t n);
void init_array(Array *a, double *data, size_t n);

MoleculeSystem init_ms(double mu, MonomerType t1, MonomerType t2, double *I1, double *I2); 
void free_ms(MoleculeSystem *ms);

int rhs(realtype t, N_Vector y, N_Vector ydot, void *data);

#endif // HAWAII_H_

#ifdef HAWAII_IMPLEMENTATION

Array create_array(size_t n) {
    Array a = {0};
    a.data = (double*) malloc(n * sizeof(double));
    assert(a.data != NULL);
    a.n = n; 
    return a;
}

void init_array(Array *a, double *data, size_t n) {
    assert(a->data != NULL);
    assert(a->n >= n);
    memcpy(a->data, data, n * sizeof(double));
}

void free_array(Array *a) {
    free(a->data);
}

MoleculeSystem init_ms(double mu, MonomerType t1, MonomerType t2, double *I1, double *I2) {
    MoleculeSystem ms = {0};
    ms.mu = mu;

    ms.m1.t = t1;

    switch (t1) {
        case ATOM: break;
        case LINEAR_MOLECULE: {
          assert(I1[0] == I1[1]);
          memcpy(ms.m1.I, I1, 2*sizeof(double));
          break;
        }
        default: TODO("init_ms");
    } 
    
    ms.m1.qp = malloc(t1 * sizeof(double));

    ms.m2.t = t2;
    switch (t2) {
        case ATOM: break;
        case LINEAR_MOLECULE: {
          assert(I2[0] == I2[1]);
          memcpy(ms.m2.I, I2, 2*sizeof(double));
          break;
        }
        default: TODO("init_ms");
    } 
    ms.m2.qp = malloc(t2 * sizeof(double));
    
    ms.Q_SIZE = t1 + t2 + 3; 
    ms.QP_SIZE = 2 * ms.Q_SIZE; 

    ms.dVdq_qp = N_VNew_Serial(ms.QP_SIZE);

    ms.dVdq = malloc(ms.Q_SIZE * sizeof(double));
    ms.intermediate_q = malloc(ms.Q_SIZE * sizeof(double));
    
    memset(ms.intermolecular_qp, 0.0, 6*sizeof(double));
    
    return ms;
}

void free_ms(MoleculeSystem *ms) {
    free(ms->m1.qp); 
    free(ms->m2.qp);
    free(ms->intermediate_q);
    free(ms->dVdq);

    N_VDestroy(ms->dVdq_qp); 
}

// перед расчетом корреляционной функции делать прикидку M0/M2
// после окончания расчета выписывать оценки M0/M2 по рассчитанной корреляционной функции

void rhsMonomer(Monomer m, double *d) {
    UNUSED(d);
    if (m.t == ATOM) return;
    TODO("rhsMonomer"); 
}


void make_qp_odd(double *q, double *qp, size_t QP_SIZE) {
    for (size_t k = 0; k < QP_SIZE; k += 2) {
        qp[k + 1] = q[k / 2];
    }
} 

void extract_q(double *qp, double *q, size_t QP_SIZE) {
    for (size_t k = 0; k < QP_SIZE; k += 2) {
        q[k / 2] = qp[k]; 
    }
}


int rhs(realtype t, N_Vector y, N_Vector ydot, void *data)
{
    UNUSED(t);
    MoleculeSystem *ms = (MoleculeSystem*) data;
  
    double Phi    = NV_Ith_S(y, IPHI); UNUSED(Phi);
    double pPhi   = NV_Ith_S(y, IPPHI);
    double Theta  = NV_Ith_S(y, ITHETA);
    double pTheta = NV_Ith_S(y, IPTHETA);
    double R      = NV_Ith_S(y, IR);
    double pR     = NV_Ith_S(y, IPR);

    double R2 = R * R;
    double R3 = R2 * R;
    double sinTheta = sin(Theta);
    double cosTheta = cos(Theta);
    double sinTheta2 = sinTheta * sinTheta;
    double sinTheta3 = sinTheta2 * sinTheta;
    
    NV_Ith_S(ydot, IR)      = pR / ms->mu;
    NV_Ith_S(ydot, IPR)     = pTheta * pTheta / (ms->mu * R3) + pPhi * pPhi / (ms->mu * R3 * sinTheta2);
    NV_Ith_S(ydot, IPHI)    = pPhi / (ms->mu * R2 * sinTheta2);
    NV_Ith_S(ydot, IPPHI)   = 0.0;
    NV_Ith_S(ydot, ITHETA)  = pTheta / (ms->mu * R2);
    NV_Ith_S(ydot, IPTHETA) = pPhi * pPhi * cosTheta / (ms->mu * R2 * sinTheta3); 
    
    double rhs_monomer1[ms->m1.t];
    rhsMonomer(ms->m1, rhs_monomer1);
    for (size_t i = 0; i < ms->m1.t; ++i) {
        NV_Ith_S(ydot, i + 6) = rhs_monomer1[i];
    }

    double rhs_monomer2[ms->m2.t];
    rhsMonomer(ms->m2, rhs_monomer2);
    for (size_t i = 0; i < ms->m2.t; ++i) {
        NV_Ith_S(ydot, i + 6 + ms->m1.t) = rhs_monomer2[i];
    }

    realtype *vdata_y = N_VGetArrayPointer(y);
    extract_q(vdata_y, ms->intermediate_q, ms->Q_SIZE);
    dpes(ms->intermediate_q, ms->dVdq);
    
    realtype *vdata_dVdq_qp = N_VGetArrayPointer(ms->dVdq_qp);
    memset(vdata_dVdq_qp, 0.0, ms->QP_SIZE);
    make_qp_odd(ms->dVdq, vdata_dVdq_qp, ms->QP_SIZE);

    N_VLinearSum(1.0, ydot, -1.0, ms->dVdq_qp, ydot); 
    
    return 0;
}

double kinetic_energy(MoleculeSystem *ms) {
    double Phi = ms->intermolecular_qp[IPHI]; UNUSED(Phi);
    double pPhi = ms->intermolecular_qp[IPPHI];
    double Theta = ms->intermolecular_qp[ITHETA];
    double pTheta = ms->intermolecular_qp[IPTHETA];
    double R = ms->intermolecular_qp[IR];
    double pR = ms->intermolecular_qp[IPR];
    
    double sinTheta = sin(Theta);
   
    double KIN1 = pR * pR / (2.0 * ms->mu) + pTheta * pTheta / (2.0 * ms->mu * R * R) + pPhi * pPhi / (2.0 * ms->mu * R * R * sinTheta * sinTheta);
    double KIN2 = 0.0;
    double KIN3 = 0.0;

    switch (ms->m1.t) {
        case ATOM: break;
        case LINEAR_MOLECULE: {
          double phi1t = ms->m1.qp[IPHI]; UNUSED(phi1t);
          double pphi1t = ms->m1.qp[IPPHI];
          double theta1t = ms->m1.qp[ITHETA];
          double ptheta1t = ms->m1.qp[IPTHETA];
          
          double sin_theta1t = sin(theta1t);
    
          KIN2 = ptheta1t * ptheta1t / (2.0 * ms->m1.I[0]) + pphi1t * pphi1t / (2.0 * ms->m1.I[1] * sin_theta1t * sin_theta1t);
          break;
        }
        default: {
          TODO("kinetic_energy");
        }
    }
    
    switch (ms->m2.t) {
        case ATOM: break;
        case LINEAR_MOLECULE: { 
          double phi2t = ms->m2.qp[IPHI]; UNUSED(phi2t);
          double pphi2t = ms->m2.qp[IPPHI];
          double theta2t = ms->m2.qp[ITHETA];
          double ptheta2t = ms->m2.qp[IPTHETA];
          
          double sin_theta2t = sin(theta2t);
    
          KIN3 = ptheta2t * ptheta2t / (2.0 * ms->m2.I[0]) + pphi2t * pphi2t / (2.0 * ms->m2.I[1] * sin_theta2t * sin_theta2t);
          break;
        }
        default: {
          TODO("kinetic_energy");
        }
    }

    return KIN1 + KIN2 + KIN3;
}

void fill_qp(MoleculeSystem *ms, Array qp)
// NOTE:
// PHI PPHI THETA PTHETA R PR 
// for monomers in the same order
{
    assert(qp.n == 6 + ms->m1.t + ms->m2.t);

    memcpy(ms->intermolecular_qp, qp.data, 6*sizeof(double));
    memcpy(ms->m1.qp, qp.data + 6, ms->m1.t*sizeof(double));
    memcpy(ms->m2.qp, qp.data + 6 + ms->m1.t, ms->m2.t*sizeof(double));
}

void print_array(Array a) {
    printf("Array[%zu] = ", a.n); 
    for (size_t i = 0; i < a.n; ++i) {
        printf("%.3f ", a.data[i]);
    }
    printf("\n"); 
}

double Hamiltonian(MoleculeSystem *ms) {
    extract_q(ms->intermolecular_qp, ms->intermediate_q, 6);
    extract_q(ms->m1.qp, ms->intermediate_q + 6/2, ms->m1.t);
    extract_q(ms->m1.qp, ms->intermediate_q + 6/2 + ms->m1.t/2, ms->m2.t);

    double V = pes(ms->intermediate_q);
    double K = kinetic_energy(ms); 

    return K + V; 
}

#endif // HAWAII_IMPLEMENTATION


