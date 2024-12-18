#include "hawaii.h"

dipolePtr dipole = NULL;

MoleculeSystem init_ms(double mu, MonomerType t1, MonomerType t2, double *I1, double *I2, size_t seed) 
{
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
    
    ms.Q_SIZE = ms.QP_SIZE / 2;
    ms.QP_SIZE = t1 + t2 + 6; 

    ms.dVdq_qp = N_VNew_Serial(ms.QP_SIZE);

    ms.dVdq = malloc(ms.Q_SIZE * sizeof(double));
    ms.intermediate_q = malloc(ms.Q_SIZE * sizeof(double));
    
    memset(ms.intermolecular_qp, 0.0, 6*sizeof(double));

    mt_seed32(seed);

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

void extract_q_and_write_into_ms(MoleculeSystem *ms) {
    extract_q(ms->intermolecular_qp, ms->intermediate_q, 6);
    extract_q(ms->m1.qp, ms->intermediate_q + 6/2, ms->m1.t);
    extract_q(ms->m1.qp, ms->intermediate_q + 6/2 + ms->m1.t/2, ms->m2.t);
}

double Hamiltonian(MoleculeSystem *ms) {
    extract_q_and_write_into_ms(ms);

    double V = pes(ms->intermediate_q);
    double K = kinetic_energy(ms); 

    // printf("R = %.5f V = %.5f\n", ms->intermolecular_qp[IR], V * HTOCM);

    return K + V; 
}

double generate_normal(double sigma) 
/*
 * Generate normally distributed variable using Box-Muller method
 */
{
    double U = mt_drand();
    double V = mt_drand();
    return sigma * sqrt(-2 * log(U)) * cos(2.0 * M_PI * V);
}

void q_generator(MoleculeSystem *ms, CalcParams *params) 
{
    double RMIN = params->sampler_Rmin;
    double RMAX = params->sampler_Rmax;
    double R3 = mt_drand()*(RMAX*RMAX*RMAX - RMIN*RMIN*RMIN) + RMIN*RMIN*RMIN;
    ms->intermolecular_qp[IR]     = pow(R3, 1.0/3.0);
    ms->intermolecular_qp[IPHI]   = mt_drand() * 2.0 * M_PI;
    ms->intermolecular_qp[ITHETA] = acos(2.0*mt_drand() - 1.0);

    switch (ms->m1.t) {
        case ATOM: break;
        case LINEAR_MOLECULE: {
          ms->m1.qp[IPHI]   = mt_drand() * 2.0 * M_PI;
          ms->m1.qp[ITHETA] = acos(2.0*mt_drand() - 1.0);
          break;
        }
        default: TODO("q_generator");
    } 
    
    switch (ms->m2.t) {
        case ATOM: break;
        case LINEAR_MOLECULE: {
          ms->m2.qp[IPHI]   = mt_drand() * 2.0 * M_PI;
          ms->m2.qp[ITHETA] = acos(2.0*mt_drand() - 1.0);
          break;
        }
        default: TODO("q_generator");
    }
}

static void p_generator_linear_molecule(Monomer *m, double Temperature)
{
    double sqrt_IIKT = sqrt(m->I[0] / HkT * Temperature);
    double sin_theta = sin(m->qp[ITHETA]); 

    double x0 = generate_normal(1.0);
    double x1 = generate_normal(1.0);

    assert(m->I[0] == m->I[1]);
    m->qp[IPPHI]   = sqrt_IIKT * sin_theta * x1;
    m->qp[IPTHETA] = sqrt_IIKT * x0;
}

static void p_generator_rotor(Monomer *m, double Temperature)
{
    double theta = m->qp[ITHETA]; 
    double psi = m->qp[IPSI]; 
    
    double sqrt_II1x = sqrt(m->I[0] / HkT * Temperature);
    double sqrt_II1y = sqrt(m->I[1] / HkT * Temperature);
    double sqrt_II1z = sqrt(m->I[2] / HkT * Temperature);

    double sin_theta, cos_theta;
    double sin_psi, cos_psi;

    sincos(theta, &sin_theta, &cos_theta);
    sincos(psi, &sin_psi, &cos_psi);
    
    double x0 = generate_normal(1.0);
    double x1 = generate_normal(1.0);
    double x2 = generate_normal(1.0);

    m->qp[IPPHI]   = sqrt_II1x * sin_theta * sin_psi * x0 + sqrt_II1y * sin_theta * cos_psi * x1 + sqrt_II1z * cos_theta * x2;
    m->qp[IPTHETA] = sqrt_II1x * cos_psi * x0 - sqrt_II1y * sin_psi * x1;
    m->qp[IPPSI]   = sqrt_II1z * x2;
}


void p_generator(MoleculeSystem *ms, double Temperature)
{
    double sqrt_MUKT = sqrt(ms->mu / HkT * Temperature);
    double sinTheta = sin(ms->intermolecular_qp[ITHETA]);

    for (size_t i = 0; i < 3; ++i) { 
        // NOTE: only 3 x-values are generated; the x-values for monomers are generated in the corresponding p_generator_... functions 
        ms->intermediate_q[i] = generate_normal(1.0);
    }

    ms->intermolecular_qp[IPPHI]   = sqrt_MUKT * ms->intermolecular_qp[IR] * sinTheta * ms->intermediate_q[0];
    ms->intermolecular_qp[IPTHETA] = sqrt_MUKT * ms->intermolecular_qp[IR]            * ms->intermediate_q[1];
    ms->intermolecular_qp[IPR]     = sqrt_MUKT * ms->intermediate_q[2];
    
    switch (ms->m1.t) {
        case ATOM: break;
        case LINEAR_MOLECULE: {
          p_generator_linear_molecule(&ms->m1, Temperature);
          break;
        }
        case ROTOR: {
          p_generator_rotor(&ms->m1, Temperature);
          break;
        }
        default: TODO("q_generator");
    }    
}

bool reject(MoleculeSystem *ms, double Temperature, double pesmin)
{
    static double PRECAUTION_FACTOR = 1.5;

    double u = mt_drand(); 

    double Mval = exp(-pesmin * HkT / Temperature); 
    double M = PRECAUTION_FACTOR * Mval;

    double proposal = exp(-kinetic_energy(ms) * HkT / Temperature); 
    double desired  = exp(-Hamiltonian(ms) * HkT / Temperature);

    double weight = desired / proposal / M;
    return weight < u;
}

void calculate_M0(MoleculeSystem *ms, CalcParams *params, double Temperature, double *m, double *q)
// Running mean/variance formulates taken from GSL 1.15
// https://github.com/ampl/gsl/blob/master/monte/plain.c 
{
    size_t counter = 0;
    size_t desired_dist = 0;
    size_t integral_counter = 0;

    double d[3];
    double fval;

    *m = 0.0;
    *q = 0.0;

    size_t print_every_nth_iteration = 1;
    switch (params->ps) {
        case FREE_AND_METASTABLE: print_every_nth_iteration = 1000000; break;
        case BOUND: print_every_nth_iteration = 1000; break;
    } 
    

    while (integral_counter < params->initialM0_npoints) {
        q_generator(ms, params);
        p_generator(ms, Temperature);

        double energy = Hamiltonian(ms);
        ++counter;

        if (!reject(ms, Temperature, params->pesmin)) {
            ++desired_dist;

            if (params->ps == FREE_AND_METASTABLE) {
                if (energy < 0.0) continue; 
            }

            if (params->ps == BOUND) {
                if (energy > 0.0) continue;
            }
    
            extract_q_and_write_into_ms(ms);
            (*dipole)(ms->intermediate_q, d);
            
            // TODO: use running mean
            fval = d[0]*d[0] + d[1]*d[1] + d[2]*d[2]; 
            double diff = fval - *m;
            *m += diff / (integral_counter + 1.0);
            *q += diff * diff * (integral_counter / (integral_counter + 1.0));
            integral_counter++;

            if (integral_counter % print_every_nth_iteration == 0) {
                printf("[calculate_M0] accumulated %zu points\n", integral_counter);
            }
        }
    } 
    
    *m = *m * ZeroCoeff * params->partial_partition_function_ratio;
    *q = sqrt(*q / integral_counter / (integral_counter - 1)) * ZeroCoeff * params->partial_partition_function_ratio;
}





