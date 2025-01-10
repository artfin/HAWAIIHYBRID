#include "hawaii.h"

#define HISTOGRAM_MAX_TPS 50

dipolePtr dipole = NULL;

MoleculeSystem init_ms(double mu, MonomerType t1, MonomerType t2, double *I1, double *I2, size_t seed) 
{
    INIT_WRANK;

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
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: {
          assert(I1[0] == I1[1]);
          memcpy(ms.m1.I, I1, 2*sizeof(double));
          break;
        }
        case ROTOR: {
          TODO("init_ms");
        } 
        case ROTOR_REQUANTIZED_ROTATION: {
          TODO("init_ms");
        }
    } 
    
    ms.m1.qp = malloc((t1 % MODULO_BASE) * sizeof(double));
    ms.m1.dVdq = malloc((t1 % MODULO_BASE)/2 * sizeof(double));

    ms.m2.t = t2;
    switch (t2) {
        case ATOM: break;
        case LINEAR_MOLECULE: {
          assert(I2[0] == I2[1]);
          memcpy(ms.m2.I, I2, 2*sizeof(double));
          break;
        }
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: {
          assert(I2[0] == I2[1]);
          memcpy(ms.m2.I, I2, 2*sizeof(double));
          break;
        }
        case ROTOR: {
          TODO("init_ms");
        } 
        case ROTOR_REQUANTIZED_ROTATION: {
          TODO("init_ms");
        }
    }

    ms.m2.qp = malloc((t2 % MODULO_BASE) * sizeof(double));
    ms.m2.dVdq = malloc((t2 % MODULO_BASE)/2 * sizeof(double));
    
    ms.QP_SIZE = (t1 % MODULO_BASE) + (t2 % MODULO_BASE) + 6; 
    ms.Q_SIZE = ms.QP_SIZE / 2;

    ms.dVdq = malloc(ms.Q_SIZE * sizeof(double));
    ms.intermediate_q = malloc(ms.Q_SIZE * sizeof(double));
    
    memset(ms.intermolecular_qp, 0.0, 6 * sizeof(double));

    if (seed > 0) { 
        mt_seed32(seed);
    } else {
        seed = mt_goodseed();
        mt_seed32(seed);
    }
     
    PRINT0("-------------------------------------------------------------------\n");
    PRINT0("    INITIALIZING MOLECULE SYSTEM %s-%s\n", monomer_type_name(t1), monomer_type_name(t2));
    PRINT0("Reduced mass of the molecule system: %.6e\n", mu);

    PRINT0("1st monomer inertia tensor [%s]: ", monomer_type_name(t1));
    switch (t1) {
        case ATOM:                                 PRINT0("\n"); break;
        case LINEAR_MOLECULE:                      PRINT0("%.3e %.3e\n", ms.m1.I[0], ms.m1.I[1]); break;
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: PRINT0("%.3e %.3e\n", ms.m1.I[0], ms.m1.I[1]); break;
        case ROTOR:                                PRINT0("%.3e %.3e %.3e\n", ms.m1.I[0], ms.m1.I[1], ms.m1.I[2]); break;
        case ROTOR_REQUANTIZED_ROTATION:           PRINT0("%.3e %.3e %.3e\n", ms.m1.I[0], ms.m1.I[1], ms.m1.I[2]); break;
    }
    
    PRINT0("2nd monomer inertia tensor [%s]: ", monomer_type_name(t2));
    switch (t2) {
        case ATOM:                                 PRINT0("\n"); break;
        case LINEAR_MOLECULE:                      PRINT0("%.3e %.3e\n", ms.m2.I[0], ms.m2.I[1]); break;
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: PRINT0("%.3e %.3e\n", ms.m2.I[0], ms.m2.I[1]); break;
        case ROTOR:                                PRINT0("%.3e %.3e %.3e\n", ms.m2.I[0], ms.m2.I[1], ms.m2.I[2]); break;
        case ROTOR_REQUANTIZED_ROTATION:           PRINT0("%.3e %.3e %.3e\n", ms.m2.I[0], ms.m2.I[1], ms.m2.I[2]); break;
    }

    PRINT0("Length of Q vector:  3 + %d + %d = %zu\n", (t1 % MODULO_BASE)/2, (t2 % MODULO_BASE)/2, ms.Q_SIZE); 
    PRINT0("Length of QP vector: 6 + %d + %d = %zu\n", (t1 % MODULO_BASE), (t2 % MODULO_BASE), ms.QP_SIZE);
    PRINT0("Generator seed is set to %zu\n", seed);
    PRINT0("-------------------------------------------------------------------\n");

    return ms;
}

void free_ms(MoleculeSystem *ms) {
    free(ms->m1.qp); 
    free(ms->m2.qp);
    free(ms->m1.dVdq);
    free(ms->m2.dVdq);

    free(ms->intermediate_q);
    free(ms->dVdq);
}

const char* monomer_type_name(MonomerType t) {
    switch (t) {
        case ATOM:                                 return "ATOM";
        case LINEAR_MOLECULE:                      return "LINEAR MOLECULE";
        case ROTOR:                                return "ROTOR";
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: return "LINEAR_MOLECULE_REQUANTIZED_ROTATION";
        case ROTOR_REQUANTIZED_ROTATION:           return "ROTOR_REQUANTIZED_ROTATION";
    }

    UNREACHABLE("monomer_type_name");
}

const char* pair_state_name(PairState ps) {
    switch (ps) {
        case FREE_AND_METASTABLE: return "FREE_AND_METASTABLE";
        case BOUND:               return "BOUND";
    }

    UNREACHABLE("pair_state_name");
}

// перед расчетом корреляционной функции делать прикидку M0/M2
// после окончания расчета выписывать оценки M0/M2 по рассчитанной корреляционной функции

void make_qp_odd(double *q, double *qp, size_t QP_SIZE) {
    for (size_t k = 0; k < QP_SIZE; k += 2) {
        qp[k + 1] = q[k / 2];
    }
} 

double find_closest_half_integer(double j) 
/*
 * requantization to:
 *    0.0, 1.5, 2.5, 3.5 ...
 * NOTE: probably should return 1.5 if j > 0.75, however it doesn't..
 */
{
    double r = 0.0;

    // to avoid the requantization to j = 0.5
    if (j < 1.0) {
        return 0.0;
    }

    while (j - r >= 1.0) {
        r += 1.0;
    }

    return r + 0.5;
}

double j_monomer(Monomer m) {
    switch (m.t) {
        case ATOM: return 0.0;
        case LINEAR_MOLECULE:
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: {
            double phi    = m.qp[IPHI]; 
            double pPhi   = m.qp[IPPHI];
            double theta  = m.qp[ITHETA];
            double pTheta = m.qp[IPTHETA];

            double jx = -pTheta * sin(phi) - pPhi * cos(phi) / tan(theta);
            double jy = pTheta * cos(phi) - pPhi * sin(phi) / tan(theta);
            double jz = pPhi;

            return sqrt(jx*jx + jy*jy + jz*jz);
        }
        case ROTOR_REQUANTIZED_ROTATION: 
        case ROTOR: {
            TODO("j_monomer");
        }
    }

    UNREACHABLE("j_monomer");
}

double torque_monomer(Monomer m)
// torque: T = dJ/dt 
{
    switch (m.t) {
        case ATOM: return 0.0;
        case LINEAR_MOLECULE:
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: {
            double phi   = m.qp[IPHI];
            double theta = m.qp[ITHETA];

            double dVdphi   = m.dVdq[IPHI / 2];
            double dVdtheta = m.dVdq[ITHETA / 2];

            double torquex = sin(phi) * dVdtheta + cos(phi) / tan(theta) * dVdphi;
            double torquey = -cos(phi) * dVdtheta + sin(phi) / tan(theta) * dVdphi;
            double torquez = -dVdphi;

            return sqrt(torquex*torquex + torquey*torquey + torquez*torquez);
        }
        case ROTOR: 
        case ROTOR_REQUANTIZED_ROTATION: {
           TODO("torque_monomer");
        }
    }

    UNREACHABLE("torque_monomer");
}

void rhsMonomer(Monomer *m, double *deriv) {
    switch (m->t) {
        case ATOM: break;
        case LINEAR_MOLECULE: {
           double pPhi   = m->qp[IPPHI];
           double Theta  = m->qp[ITHETA];
           double pTheta = m->qp[IPTHETA];

           double sin_theta = sin(Theta);
           double cos_theta = cos(Theta);

           deriv[IPHI]    = pPhi / m->I[0] / sin_theta / sin_theta;
           deriv[IPPHI]   = -m->dVdq[IPHI/2];
           deriv[ITHETA]  = pTheta / m->I[0]; 
           deriv[IPTHETA] = pPhi * pPhi * cos_theta / m->I[0] / sin_theta / sin_theta / sin_theta - m->dVdq[ITHETA/2]; 
           
           break;                                                    
        }
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: {
           if (m->apply_requantization) {
                double j    = j_monomer(*m);
                double jreq = find_closest_half_integer(j);

                double scaling_factor = 0.0;
                if (j > 1e-15) {
                    scaling_factor = jreq / j; 
                }

                m->qp[IPPHI]   *= scaling_factor;
                m->qp[IPTHETA] *= scaling_factor;
           }

           double pPhi   = m->qp[IPPHI];
           double Theta  = m->qp[ITHETA];
           double pTheta = m->qp[IPTHETA];

           double sin_theta = sin(Theta);
           double cos_theta = cos(Theta);

           deriv[IPHI]    = pPhi / m->I[0] / sin_theta / sin_theta;
           deriv[IPPHI]   = -m->dVdq[IPHI/2];
           deriv[ITHETA]  = pTheta / m->I[0]; 
           deriv[IPTHETA] = pPhi * pPhi * cos_theta / m->I[0] / sin_theta / sin_theta / sin_theta - m->dVdq[ITHETA/2]; 
    
           break;                                                    
        }
        case ROTOR: {
          TODO("rhsMonomer");
        }
        case ROTOR_REQUANTIZED_ROTATION: {
          TODO("rhsMonomer");
        }
    } 
}

void extract_dVdq_and_write_into_monomers(MoleculeSystem *ms) {
    memcpy(ms->m1.dVdq, ms->dVdq + 3,                              (ms->m1.t % MODULO_BASE)/2 * sizeof(double));
    memcpy(ms->m2.dVdq, ms->dVdq + 3 + (ms->m1.t % MODULO_BASE)/2, (ms->m2.t % MODULO_BASE)/2 * sizeof(double));
}

int rhs(realtype t, N_Vector y, N_Vector ydot, void *data)
{
    UNUSED(t);

    assert(data != NULL); 
    MoleculeSystem *ms = (MoleculeSystem*) data;
    put_qp_into_ms(ms, (Array){.data = N_VGetArrayPointer(y), .n = ms->QP_SIZE});
    
    extract_q_and_write_into_ms(ms);
    dpes(ms->intermediate_q, ms->dVdq);
    extract_dVdq_and_write_into_monomers(ms);
    
    double Phi    = ms->intermolecular_qp[IPHI]; UNUSED(Phi);
    double pPhi   = ms->intermolecular_qp[IPPHI];
    double Theta  = ms->intermolecular_qp[ITHETA];
    double pTheta = ms->intermolecular_qp[IPTHETA];
    double R      = ms->intermolecular_qp[IR];
    double pR     = ms->intermolecular_qp[IPR];

    double R2 = R * R;
    double R3 = R2 * R;
    double sinTheta = sin(Theta);
    double cosTheta = cos(Theta);
    double sinTheta2 = sinTheta * sinTheta;
    double sinTheta3 = sinTheta2 * sinTheta;
    
    NV_Ith_S(ydot, IR)      = pR / ms->mu;
    NV_Ith_S(ydot, IPR)     = pTheta * pTheta / (ms->mu * R3) + pPhi * pPhi / (ms->mu * R3 * sinTheta2) - ms->dVdq[IR/2];
    NV_Ith_S(ydot, IPHI)    = pPhi / (ms->mu * R2 * sinTheta2);
    NV_Ith_S(ydot, IPPHI)   = -ms->dVdq[IPHI/2];
    NV_Ith_S(ydot, ITHETA)  = pTheta / (ms->mu * R2);
    NV_Ith_S(ydot, IPTHETA) = pPhi * pPhi * cosTheta / (ms->mu * R2 * sinTheta3) - ms->dVdq[ITHETA/2]; 
    
    double rhs_monomer1[ms->m1.t % MODULO_BASE];
    rhsMonomer(&ms->m1, rhs_monomer1);
    for (size_t i = 0; i < ms->m1.t % MODULO_BASE; ++i) {
        NV_Ith_S(ydot, i + 6) = rhs_monomer1[i];
    }

    double rhs_monomer2[ms->m2.t % MODULO_BASE];
    rhsMonomer(&ms->m2, rhs_monomer2);
    for (size_t i = 0; i < ms->m2.t % MODULO_BASE; ++i) {
        NV_Ith_S(ydot, i + 6 + (ms->m1.t % MODULO_BASE)) = rhs_monomer2[i];
    }
    
    //for (size_t i = 0; i < 10; ++i) {
    //    printf("y(%zu) = %.15e\n", i, NV_Ith_S(y, i)); 
    //}

    //for (size_t i = 0; i < 10; ++i) {
    //    printf("ydot(%zu) = %.15e\n", i, NV_Ith_S(ydot, i));
    //}
    //assert(false);

    return 0;
}

Array compute_numerical_rhs(MoleculeSystem *ms) 
{
    Array derivatives = create_array(ms->QP_SIZE);
    Array qp = create_array(ms->QP_SIZE);
   
    double step = 1e-7;

    for (size_t i = 0; i < ms->QP_SIZE; ++i) {
        get_qp_from_ms(ms, &qp);

        if (i % 2 == 0) {
            double c = qp.data[i + 1];

            qp.data[i + 1] = c + step;
            put_qp_into_ms(ms, qp);
            double Ep = Hamiltonian(ms);
            
            qp.data[i + 1] = c - step;
            put_qp_into_ms(ms, qp);
            double Em = Hamiltonian(ms);
            
            derivatives.data[i] = (Ep - Em)/(2.0*step);
        } else {
            double c = qp.data[i - 1];

            qp.data[i - 1] = c + step;
            put_qp_into_ms(ms, qp);
            double Ep = Hamiltonian(ms);

            qp.data[i - 1] = c - step;
            put_qp_into_ms(ms, qp);
            double Em = Hamiltonian(ms);

            derivatives.data[i] = -(Ep - Em) / (2.0*step);
        }
    }

    free_array(&qp);

    return derivatives;
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
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: 
        case LINEAR_MOLECULE: {
          double phi1t = ms->m1.qp[IPHI]; UNUSED(phi1t);
          double pphi1t = ms->m1.qp[IPPHI];
          double theta1t = ms->m1.qp[ITHETA];
          double ptheta1t = ms->m1.qp[IPTHETA];
          
          double sin_theta1t = sin(theta1t);
           
          KIN2 = ptheta1t * ptheta1t / (2.0 * ms->m1.I[0]) + pphi1t * pphi1t / (2.0 * ms->m1.I[1] * sin_theta1t * sin_theta1t);
          break;
        }
        case ROTOR_REQUANTIZED_ROTATION: 
        case ROTOR: { 
          TODO("kinetic_energy");
        }
    }
    
    switch (ms->m2.t) {
        case ATOM: break;
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: 
        case LINEAR_MOLECULE: { 
          double phi2t = ms->m2.qp[IPHI]; UNUSED(phi2t);
          double pphi2t = ms->m2.qp[IPPHI];
          double theta2t = ms->m2.qp[ITHETA];
          double ptheta2t = ms->m2.qp[IPTHETA];
          
          double sin_theta2t = sin(theta2t);
    
          KIN3 = ptheta2t * ptheta2t / (2.0 * ms->m2.I[0]) + pphi2t * pphi2t / (2.0 * ms->m2.I[1] * sin_theta2t * sin_theta2t);
          break;
        }
        case ROTOR_REQUANTIZED_ROTATION: 
        case ROTOR: { 
          TODO("kinetic_energy");
        }
    }


    return KIN1 + KIN2 + KIN3;
}

void put_qp_into_ms(MoleculeSystem *ms, Array qp)
// NOTE:
// PHI PPHI THETA PTHETA R PR 
// for monomers in the same order
{
    assert(qp.n == 6 + (ms->m1.t%MODULO_BASE) + (ms->m2.t%MODULO_BASE));

    memcpy(ms->intermolecular_qp, qp.data,                          6*sizeof(double));
    memcpy(ms->m1.qp,             qp.data+6,                        (ms->m1.t%MODULO_BASE)*sizeof(double));
    memcpy(ms->m2.qp,             qp.data+6+(ms->m1.t%MODULO_BASE), (ms->m2.t%MODULO_BASE)*sizeof(double));
}

void get_qp_from_ms(MoleculeSystem *ms, Array *qp) {
    assert(qp->n == ms->QP_SIZE);

    memcpy(qp->data,                          ms->intermolecular_qp, 6*sizeof(double));
    memcpy(qp->data+6,                        ms->m1.qp,             (ms->m1.t%MODULO_BASE)*sizeof(double));
    memcpy(qp->data+6+(ms->m1.t%MODULO_BASE), ms->m2.qp,             (ms->m2.t%MODULO_BASE)*sizeof(double));
}

void extract_q(double *qp, double *q, size_t QP_SIZE) {
    for (size_t k = 0; k < QP_SIZE; k += 2) {
        q[k / 2] = qp[k]; 
    }
}

void extract_q_and_write_into_ms(MoleculeSystem *ms) {
    extract_q(ms->intermolecular_qp, ms->intermediate_q,                              6);
    extract_q(ms->m1.qp,             ms->intermediate_q+6/2,                          ms->m1.t%MODULO_BASE);
    extract_q(ms->m1.qp,             ms->intermediate_q+6/2+(ms->m1.t%MODULO_BASE)/2, ms->m2.t%MODULO_BASE);
}

void invert_momenta(MoleculeSystem *ms) {
    for (size_t i = 1; i < 6; i += 2) {
        ms->intermolecular_qp[i] = -ms->intermolecular_qp[i];
    }

    for (size_t i = 1; i < (ms->m1.t%MODULO_BASE); ++i) {
        ms->m1.qp[i] = -ms->m1.qp[i];
    }
    
    for (size_t i = 1; i < (ms->m2.t%MODULO_BASE); ++i) {
        ms->m2.qp[i] = -ms->m2.qp[i];
    }
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
    assert(params->sampler_Rmin > 0);
    assert(params->sampler_Rmax > 0);

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

    // NOTE: only 3 x-values are generated; the x-values for monomers are generated in the corresponding p_generator_... functions 
    for (size_t i = 0; i < 3; ++i) { 
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
// Running mean/variance formulas taken from GSL 1.15
// https://github.com/ampl/gsl/blob/master/monte/plain.c 
{
    assert(params->initialM0_npoints > 0);
    assert(fabs(params->pesmin) > 1e-15);

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
        case BOUND:               print_every_nth_iteration = 1000;    break;
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
                double M0_est = *m * ZeroCoeff * params->partial_partition_function_ratio;
                double M0std_est = sqrt(*q / integral_counter / (integral_counter - 1)) * ZeroCoeff * params->partial_partition_function_ratio;
                printf("[calculate_M0] %zu/%zu points: \t M0 = %.5e +/- %.5e\n", integral_counter, params->initialM0_npoints, M0_est, M0std_est);
            }
        }
    } 
    
    *m = *m * ZeroCoeff * params->partial_partition_function_ratio;
    *q = sqrt(*q / integral_counter / (integral_counter - 1)) * ZeroCoeff * params->partial_partition_function_ratio;
}

#ifdef USE_MPI
void mpi_calculate_M0(MoleculeSystem *ms, CalcParams *params, double Temperature, double *m, double *q)
{
    assert(params->initialM0_npoints > 0);
    assert(fabs(params->pesmin) > 1e-15);

    INIT_WRANK;
    INIT_WSIZE;

    size_t counter = 0;
    size_t desired_dist = 0;
    size_t integral_counter = 0;

    double d[3];
    double fval;

    size_t print_every_nth_iteration = 1;
    switch (params->ps) {
        case FREE_AND_METASTABLE: print_every_nth_iteration = 1000000; break;
        case BOUND:               print_every_nth_iteration = 1000;    break;
    } 
    
    size_t local_npoints = params->initialM0_npoints / _wsize;
    
    *m = 0.0;
    *q = 0.0;
    double ml = 0.0, ql = 0.0;

    while (integral_counter < local_npoints) {
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
            double diff = fval - ml;
            ml += diff / (integral_counter + 1.0);
            ql += diff * diff * (integral_counter / (integral_counter + 1.0));
            integral_counter++;

            if (integral_counter % print_every_nth_iteration == 0) {
                double M0_est    = ml * ZeroCoeff * params->partial_partition_function_ratio;
                double M0std_est = sqrt(ql / integral_counter / (integral_counter - 1)) * ZeroCoeff * params->partial_partition_function_ratio;
                PRINT0("[mpi_calculate_M0] %zu/%zu points: \t M0 = %.5e +/- %.5e\n", _wsize*integral_counter, params->initialM0_npoints, M0_est, M0std_est);
            }
        }
    } 

    MPI_Allreduce(MPI_IN_PLACE, &ml, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &ql, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
     
    *m = (ml/_wsize) * ZeroCoeff * params->partial_partition_function_ratio;
    *q = sqrt((ql/_wsize) / integral_counter / (integral_counter - 1)) * ZeroCoeff * params->partial_partition_function_ratio;
}
#endif // USE_MPI


#include "trajectory.h"

void track_turning_points(Tracker *tr, double R)
{
    tr->before2 = tr->before;
    tr->before  = tr->current;
    tr->current = R; 

    if (!tr->empty) {
        // local maxima
        if ((tr->before > tr->before2) && (tr->before > tr->current)) {
            tr->turning_points++;
        }

        // local minima
        if ( (tr->before < tr->before2) && (tr->before < tr->current)) { 
            tr->turning_points++;
        }
    }

    tr->called++;
    if (tr->called > 1) tr->empty = false; 
} 

int correlation_eval(MoleculeSystem *ms, Trajectory *traj, CalcParams *params, double *crln, size_t *tps)
// TODO: Use temporary arena instead of malloc 
{
    double *correlation_forw = malloc(params->MaxTrajectoryLength * sizeof(double));
    memset(correlation_forw, 0.0, params->MaxTrajectoryLength * sizeof(double));

    double *correlation_back = malloc(params->MaxTrajectoryLength * sizeof(double));
    memset(correlation_back, 0.0, params->MaxTrajectoryLength * sizeof(double));
    
    memset(crln, 0, params->MaxTrajectoryLength * sizeof(double));
            
    double dip0[3], dipt[3];
    extract_q_and_write_into_ms(ms);
    (*dipole)(ms->intermediate_q, dip0);
   
    correlation_forw[0] = dip0[0]*dip0[0] + dip0[1]*dip0[1] + dip0[2]*dip0[2]; 
    correlation_back[0] = dip0[0]*dip0[0] + dip0[1]*dip0[1] + dip0[2]*dip0[2]; 
    
    Array qp = create_array(ms->QP_SIZE);
    get_qp_from_ms(ms, &qp);
    set_initial_condition(traj, qp);
    
    int status = 0;
    
    double t = 0.0;
    double tout = params->sampling_time;
   
    Tracker tr = {
      .before2 = qp.data[IR],
      .before  = qp.data[IR],
      .current = qp.data[IR],
      .empty   = true,
    };

    /*
     * We start step_counter from 1 so that correlation value after the first integration step
     * will go into correlation_forw[1] 
     */
    for (size_t step_counter = 1; step_counter < params->MaxTrajectoryLength; ++step_counter, tout += params->sampling_time)
    {
        status = make_step(traj, tout, &t);
        if (status) {
            fprintf(stderr, "CVODE ERROR: status = %d\n", status);
            break;
        }

        put_qp_into_ms(ms, (Array){.data = N_VGetArrayPointer(traj->y), .n = ms->QP_SIZE});
        extract_q_and_write_into_ms(ms);
        (*dipole)(ms->intermediate_q, dipt);

        correlation_forw[step_counter] = dip0[0]*dipt[0] + dip0[1]*dipt[1] + dip0[2]*dipt[2]; 

        track_turning_points(&tr, ms->intermolecular_qp[IR]);

        if (ms->intermolecular_qp[IR] > params->Rcut) break;
    }
    
    put_qp_into_ms(ms, qp);
    invert_momenta(ms);
    get_qp_from_ms(ms, &qp);
    set_initial_condition(traj, qp); // re-initialization of the CVode happens here 
   
    
    t = 0.0;
    tout = params->sampling_time;
   
    tr.before2 = qp.data[IR];
    tr.before  = qp.data[IR];
    tr.current = qp.data[IR];
    tr.called  = 0;
    tr.empty   = true;

    for (size_t step_counter = 1; step_counter < params->MaxTrajectoryLength; ++step_counter, tout += params->sampling_time)
    {
        status = make_step(traj, tout, &t);
        if (status) {
            fprintf(stderr, "CVODE ERROR: status = %d\n", status);
            break;
        }

        put_qp_into_ms(ms, (Array){.data = N_VGetArrayPointer(traj->y), .n = ms->QP_SIZE});
        extract_q_and_write_into_ms(ms);
        (*dipole)(ms->intermediate_q, dipt);

        correlation_back[step_counter] = dip0[0]*dipt[0] + dip0[1]*dipt[1] + dip0[2]*dipt[2]; 
        
        track_turning_points(&tr, ms->intermolecular_qp[IR]);

        if (ms->intermolecular_qp[IR] > params->Rcut) break;
    }

    for (size_t i = 0; i < params->MaxTrajectoryLength; ++i) {
        crln[i] = 0.5 * (correlation_forw[i] + correlation_back[i]);
    } 

    *tps = tr.turning_points;

    free_array(&qp);
   
    free(correlation_forw);
    free(correlation_back); 

    return status;
}

#ifdef USE_MPI
CFnc calculate_correlation_and_save(MoleculeSystem *ms, CalcParams *params, double Temperature)
{
    assert(dipole != NULL);

    assert(params->MaxTrajectoryLength > 0);
    assert(params->Rcut > 0);
    assert(params->sampling_time > 0);
    assert(params->total_trajectories > 0);
    assert(params->partial_partition_function_ratio > 0);
    assert(params->cvode_tolerance > 0);
    assert(params->niterations >= 1);
    assert(params->cf_filename != NULL);

    INIT_WRANK;
    INIT_WSIZE;
  
    FILE *fd = NULL; 
    if (_wrank == 0) {
        fd = fopen(params->cf_filename, "w");
        if (fd == NULL) { 
            printf("ERROR: Could not open '%s' for writing! Exiting...\n", params->cf_filename);
            exit(1);
        }
    } 

    size_t local_ntrajectories = params->total_trajectories / params->niterations / _wsize;
   
    double *crln       = malloc(params->MaxTrajectoryLength * sizeof(double));
    double *local_crln = malloc(params->MaxTrajectoryLength * sizeof(double));
    
    CFnc total_crln = {
        .t     = linspace(0.0, params->sampling_time*(params->MaxTrajectoryLength-1), params->MaxTrajectoryLength),
        .data  = malloc(params->MaxTrajectoryLength * sizeof(double)),
        .len   = params->MaxTrajectoryLength,
        .ntraj = 0,
        .T     = Temperature, 
    }; 
    memset(total_crln.data, 0, params->MaxTrajectoryLength * sizeof(double)); 
    
    CFnc total_crln_iter = {
        .t     = linspace(0.0, params->sampling_time*(params->MaxTrajectoryLength-1), params->MaxTrajectoryLength),
        .data  = malloc(params->MaxTrajectoryLength * sizeof(double)),
        .len   = params->MaxTrajectoryLength,
        .ntraj = 0,
        .T     = Temperature, 
    }; 
    memset(total_crln_iter.data, 0, params->MaxTrajectoryLength * sizeof(double)); 
    
    size_t print_every_nth_iteration = 1;
    switch (params->ps) {
        case FREE_AND_METASTABLE: print_every_nth_iteration = 1000000; break;
        case BOUND:               print_every_nth_iteration = 1000;    break;
    } 
    
    Trajectory traj = init_trajectory(ms, params->cvode_tolerance);

    gsl_histogram *tps_hist = NULL;
    if (params->ps == FREE_AND_METASTABLE) {
        size_t nbins = HISTOGRAM_MAX_TPS;
        tps_hist = gsl_histogram_alloc(nbins);
        gsl_histogram_set_ranges_uniform(tps_hist, 0, HISTOGRAM_MAX_TPS);
    }

    PRINT0("\n\n"); 
    PRINT0("------------------------------------------------------------------------\n");
    PRINT0("Calculating single correlation function at T = %.2f using following parameters:\n", Temperature);
    PRINT0("    pair state (pair_state):                                             %s\n",     pair_state_name(params->ps));
    PRINT0("    trajectories to be calculated (total_trajectories):                  %zu\n",    params->total_trajectories);
    PRINT0("    # of iterations that the calculation is divided into (niterations):  %zu\n",    params->niterations);
    PRINT0("    maximum length of trajectory (MaxTrajectoryLength):                  %zu\n",    params->MaxTrajectoryLength);
    PRINT0("    partial partition function (partial_partition_function_ratio):       %.6e\n",   params->partial_partition_function_ratio);
    PRINT0("    sampling time of dipole on trajectory (sampling_time):               %.2f\n",   params->sampling_time);
    PRINT0("    maximum intermolecular distance on trajectory (Rcut):                %.2f\n",   params->Rcut);
    PRINT0("    CVode tolerance:                                                     %.3e\n\n", params->cvode_tolerance);
    PRINT0("------------------------------------------------------------------------\n");
    PRINT0("\n\n"); 
    PRINT0("Running preliminary calculation of M0 to test the sampler and dipole function...\n");
    PRINT0("The estimate will be based on %zu points\n\n", params->initialM0_npoints); 

    double prelim_M0, prelim_M0std;
    mpi_calculate_M0(ms, params, Temperature, &prelim_M0, &prelim_M0std);
    PRINT0("M0 = %.10e +/- %.10e [%.10e ... %.10e]\n", prelim_M0, prelim_M0std, prelim_M0 - prelim_M0std, prelim_M0 + prelim_M0std);
    PRINT0("Error: %.3f%%\n", prelim_M0std/prelim_M0 * 100.0);
    
    for (size_t iter = 0; iter < params->niterations; ++iter) 
    {
        size_t counter = 0;
        size_t desired_dist = 0;
        size_t integral_counter = 0;
        size_t tps = 0; 

        while (integral_counter < local_ntrajectories) 
        {
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

                int status = correlation_eval(ms, &traj, params, crln, &tps); 
                if (status == -1) continue;

                if (params->ps == FREE_AND_METASTABLE) {
                    gsl_histogram_increment(tps_hist, tps);
                }

                for (size_t i = 0; i < params->MaxTrajectoryLength; ++i) {
                    local_crln[i] += params->partial_partition_function_ratio * crln[i];
                }

                integral_counter++;

                if (integral_counter % print_every_nth_iteration == 0) {
                    printf("[%d - calculate_correlation] accumulated %zu points\n", _wrank, integral_counter);
                }
            }
        }

        MPI_Allreduce(local_crln, total_crln_iter.data, params->MaxTrajectoryLength, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
       
        for (size_t i = 0; i < params->MaxTrajectoryLength; ++i) {
            total_crln.data[i] += total_crln_iter.data[i];
        }
        total_crln.ntraj += local_ntrajectories * _wsize;

        memset(local_crln, 0, params->MaxTrajectoryLength * sizeof(double));
        memset(total_crln_iter.data, 0, params->MaxTrajectoryLength * sizeof(double));

        PRINT0("ITERATION %zu/%zu: accumulated %zu trajectories. Saving the temporary result to '%s'\n", iter, params->niterations, total_crln.ntraj, params->cf_filename);
        double M0_crln_est = total_crln.data[0] / total_crln.ntraj * ZeroCoeff;
        PRINT0("M0 ESTIMATE FROM CF: %.5e, PRELIMINARY M0 ESTIMATE: %.5e, diff: %.3f%%\n\n", M0_crln_est, prelim_M0, (M0_crln_est - prelim_M0)/prelim_M0*100.0);

        if (_wrank == 0) save_correlation_function(fd, total_crln, params);
    }
 
    if (tps_hist != NULL) {
        gsl_histogram_free(tps_hist);
    }
    
    free(crln);
    free(local_crln);
    free_cfnc(total_crln_iter);

    return total_crln; 
}
#endif // USE_MPI

int assert_float_is_equal_to(double estimate, double true_value, double abs_tolerance) {
    INIT_WRANK;
    INIT_WSIZE;

    if ((estimate > (true_value - abs_tolerance)) && (estimate < (true_value + abs_tolerance))) {
        PRINT0("\033[32mASSERTION PASSED:\033[0m Estimate lies within expected bounds from true value!\n");
        return 0; 
    } else {
        PRINT0("\033[31mASSERTION FAILED:\033[0m\n");
        PRINT0("ERROR: Estimate lies outside expected bounds from true value!\n");
        PRINT0("Expected bounds: %.5e...%.5e and received %.5e\n", true_value - abs_tolerance, true_value + abs_tolerance, estimate);
        return 1; 
    }

    UNREACHABLE("");
}

double* linspace(double start, double end, size_t n) {
    double step = (end - start)/(n - 1);

    if (n == 1) {
        assert(start == end);
        step = end - start;
    }
    
    double *v = (double*) malloc(n * sizeof(double)); 
    for (size_t i = 0; i < n; ++i) {
        v[i] = start + step * i;
    }

    return v;
}

void save_correlation_function(FILE *fd, CFnc crln, CalcParams *params)
{
    fprintf(fd, "# HAWAII HYBRID v0.1\n");
    fprintf(fd, "# TEMPERATURE: %.2f\n", crln.T);
    fprintf(fd, "# PAIR STATE: %s\n", pair_state_name(params->ps));
    fprintf(fd, "# AVERAGE OVER %zu TRAJECTORIES\n", crln.ntraj); 
    fprintf(fd, "# MAXIMUM TRAJECTORY LENGTH: %zu\n", params->MaxTrajectoryLength);
    fprintf(fd, "# PARTIAL PARTITION FUNCTION: %.2e\n", params->partial_partition_function_ratio);
    fprintf(fd, "# CVODE TOLERANCE: %.2e\n", params->cvode_tolerance);

    for (size_t i = 0; i < crln.len; ++i) {
        fprintf(fd, "%.10f %.10e\n", crln.t[i], crln.data[i] / crln.ntraj *ALU*ALU*ALU);
    }
}

double analytic_full_partition_function_by_V(MoleculeSystem *ms, double T)
{
    double pf_analytic = 0.0;

    if ((ms->m1.t == ATOM) && (ms->m2.t == ATOM)) {
        TODO("analytic_full_partition_function_by_V");  
    } else if ((ms->m1.t == LINEAR_MOLECULE) && (ms->m2.t == ATOM)) {
        pf_analytic = 4.0 * M_PI * pow(2.0 * M_PI * T / HkT, 2.5) * pow(ms->mu, 1.5) * ms->m1.I[0]; 
    }

    return pf_analytic;
}

gsl_histogram* gsl_histogram_extend_right(gsl_histogram* h)
{
    size_t nbins = h->n;
    double add_bins = 10;
    
    double xmin = h->range[0];
    double xmax = h->range[nbins];
    double dx = h->range[1] - h->range[0];
    
    double new_xmax = xmax + add_bins*dx; 
    
    gsl_histogram *new_h = gsl_histogram_alloc(nbins + add_bins);
    gsl_histogram_set_ranges_uniform(new_h, xmin, new_xmax);
        
    size_t nc = 0; // cursor over the new histogram
    for (size_t i = 0; i < nbins; ++i) {
        new_h->bin[nc++] = gsl_histogram_get(h, i); 
    }

    gsl_histogram_free(h);

    return new_h;
}

void free_cfnc(CFnc cf) {
    free(cf.t);
    free(cf.data);
}

