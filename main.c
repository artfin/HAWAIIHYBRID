#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>

#define IR 0 
#define IPR 1    
#define IPHI 2   
#define IPPHI 3  
#define ITHETA 4 
#define IPTHETA 5

typedef enum {
    ATOM = 0, 
    LINEAR_MOLECULE = 4, 
    ROTOR = 6,
} MonomerType;

typedef struct { 
    MonomerType t;
    double *I;
    double *qp;
} Monomer; 

typedef struct {
    Monomer m1;
    Monomer m2;
    double mu;

    size_t Q_SIZE;
    size_t QP_SIZE;
    
    double *intermediate_q;
    double *dVdq;
    N_Vector dVdq_qp;
} MoleculeSystem;

MoleculeSystem init_ms(MonomerType t1, MonomerType t2, double *I1, double *I2) {
    MoleculeSystem ms = {0};

    ms.m1.t = t1;
    ms.m1.I = I1;
    ms.m1.qp = malloc(t1 * sizeof(double));

    ms.m2.t = t2;
    ms.m2.I = I2;
    ms.m2.qp = malloc(t2 * sizeof(double));
    
    ms.Q_SIZE = t1 + t2 + 3; 
    ms.QP_SIZE = 2 * ms.Q_SIZE; 

    ms.dVdq_qp = N_VNew_Serial(ms.QP_SIZE);

    ms.dVdq = malloc(ms.Q_SIZE * sizeof(double));
    ms.intermediate_q = malloc(ms.Q_SIZE * sizeof(double));
    
    return ms;
}

void free_ms(MoleculeSystem *ms) {
    free(ms->m1.qp); 
    free(ms->m2.qp);
    free(ms->intermediate_q);
    free(ms->dVdq);

    N_VDestroy(ms->dVdq_qp); 
}

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

// перед расчетом корреляционной функции делать прикидку M0/M2
// после окончания расчета выписывать оценки M0/M2 по рассчитанной корреляционной функции

void rhsMonomer(Monomer m, double *d) {
    if (m.t == ATOM) return; 
}

double pes(double *q);
double dpes(double *q, double *dVdq);

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
    MoleculeSystem *ms = (MoleculeSystem*) data;
  
    double R      = NV_Ith_S(y, IR);
    double pR     = NV_Ith_S(y, IPR);
    double Phi    = NV_Ith_S(y, IPHI);
    double pPhi   = NV_Ith_S(y, IPPHI);
    double Theta  = NV_Ith_S(y, ITHETA);
    double pTheta = NV_Ith_S(y, IPTHETA);

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


int main()
{
    double rhs[0] = {};
    printf("%p\n", rhs);
    printf("%lf\n", *rhs);

    return 0;
}
