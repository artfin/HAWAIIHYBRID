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
    
    printf("-------------------------------------------------------------------\n");
    printf("    INITIALIZING MOLECULE SYSTEM %s-%s\n", monomer_type_name(t1), monomer_type_name(t2));
    printf("Reduced mass of the molecule system: %.6e\n", mu);

    printf("1st monomer inertia tensor [%s]: ", monomer_type_name(t1));
    switch (t1) {
        case ATOM:                                 printf("\n"); break;
        case LINEAR_MOLECULE:                      printf("%.3e %.3e\n", ms.m1.I[0], ms.m1.I[1]); break;
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: printf("%.3e %.3e\n", ms.m1.I[0], ms.m1.I[1]); break;
        case ROTOR:                                printf("%.3e %.3e %.3e\n", ms.m1.I[0], ms.m1.I[1], ms.m1.I[2]); break;
        case ROTOR_REQUANTIZED_ROTATION:           printf("%.3e %.3e %.3e\n", ms.m1.I[0], ms.m1.I[1], ms.m1.I[2]); break;
    }
    
    printf("2nd monomer inertia tensor [%s]: ", monomer_type_name(t2));
    switch (t2) {
        case ATOM:                                 printf("\n"); break;
        case LINEAR_MOLECULE:                      printf("%.3e %.3e\n", ms.m2.I[0], ms.m2.I[1]); break;
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: printf("%.3e %.3e\n", ms.m2.I[0], ms.m2.I[1]); break;
        case ROTOR:                                printf("%.3e %.3e %.3e\n", ms.m2.I[0], ms.m2.I[1], ms.m2.I[2]); break;
        case ROTOR_REQUANTIZED_ROTATION:           printf("%.3e %.3e %.3e\n", ms.m2.I[0], ms.m2.I[1], ms.m2.I[2]); break;
    }

    printf("Length of Q vector:  3 + %d + %d = %zu\n", (t1 % MODULO_BASE)/2, (t2 % MODULO_BASE)/2, ms.Q_SIZE); 
    printf("Length of QP vector: 6 + %d + %d = %zu\n", (t1 % MODULO_BASE), (t2 % MODULO_BASE), ms.QP_SIZE);
    printf("Generator seed is set to %zu\n", seed);
    printf("-------------------------------------------------------------------\n");

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


// перед расчетом корреляционной функции делать прикидку M0/M2
// после окончания расчета выписывать оценки M0/M2 по рассчитанной корреляционной функции

void make_qp_odd(double *q, double *qp, size_t QP_SIZE) {
    for (size_t k = 0; k < QP_SIZE; k += 2) {
        qp[k + 1] = q[k / 2];
    }
} 


void rhsMonomer(Monomer m, double *deriv) {
    switch (m.t) {
        case ATOM: break;
        case LINEAR_MOLECULE: {
           double pPhi   = m.qp[IPPHI];
           double Theta  = m.qp[ITHETA];
           double pTheta = m.qp[IPTHETA];

           double sin_theta = sin(Theta);
           double cos_theta = cos(Theta);

           deriv[IPHI]    = pPhi / m.I[0] / sin_theta / sin_theta;
           deriv[IPPHI]   = -m.dVdq[IPHI/2];
           deriv[ITHETA]  = pTheta / m.I[0]; 
           deriv[IPTHETA] = pPhi * pPhi * cos_theta / m.I[0] / sin_theta / sin_theta / sin_theta - m.dVdq[ITHETA/2]; 
           
           break;                                                    
        }
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: {
           double pPhi   = m.qp[IPPHI];
           double Theta  = m.qp[ITHETA];
           double pTheta = m.qp[IPTHETA];

           double sin_theta = sin(Theta);
           double cos_theta = cos(Theta);

           deriv[IPHI]    = pPhi / m.I[0] / sin_theta / sin_theta;
           deriv[IPPHI]   = -m.dVdq[IPHI/2];
           deriv[ITHETA]  = pTheta / m.I[0]; 
           deriv[IPTHETA] = pPhi * pPhi * cos_theta / m.I[0] / sin_theta / sin_theta / sin_theta - m.dVdq[ITHETA/2]; 
        
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

double torque_monomer(MoleculeSystem *ms, size_t monomer_index) 
{
    assert((monomer_index == 0) || (monomer_index == 1));

    Monomer *m = NULL;
    if (monomer_index == 0) { 
        m = &ms->m1;

        switch (m->t) {
          case ATOM: return 0.0;
          case LINEAR_MOLECULE:
          case LINEAR_MOLECULE_REQUANTIZED_ROTATION: {
            double phi   = m->qp[IPHI];
            double theta = m->qp[ITHETA];
        
            extract_q_and_write_into_ms(ms);
            dpes(ms->intermediate_q, ms->dVdq);

            double dVdphi   = ms->dVdq[6/2 + IPHI];
            double dVdtheta = ms->dVdq[6/2 + ITHETA];

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
    }

    if (monomer_index == 1) {
        m = &ms->m2;

        switch (m->t) {
          case ATOM: return 0.0;
          case LINEAR_MOLECULE:
          case LINEAR_MOLECULE_REQUANTIZED_ROTATION: {
            double phi   = m->qp[IPHI];
            double theta = m->qp[ITHETA];
        
            extract_q_and_write_into_ms(ms);
            dpes(ms->intermediate_q, ms->dVdq);

            double dVdphi   = ms->dVdq[6/2 + (ms->m1.t % MODULO_BASE) + IPHI];
            double dVdtheta = ms->dVdq[6/2 + (ms->m1.t % MODULO_BASE) + ITHETA];

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
    }

    UNREACHABLE("torque_monomer");
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
    rhsMonomer(ms->m1, rhs_monomer1);
    for (size_t i = 0; i < ms->m1.t % MODULO_BASE; ++i) {
        NV_Ith_S(ydot, i + 6) = rhs_monomer1[i];
    }

    double rhs_monomer2[ms->m2.t % MODULO_BASE];
    rhsMonomer(ms->m2, rhs_monomer2);
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
    extract_q(ms->intermolecular_qp, ms->intermediate_q, 6);
    extract_q(ms->m1.qp, ms->intermediate_q + 6/2, ms->m1.t % MODULO_BASE);
    extract_q(ms->m1.qp, ms->intermediate_q + 6/2 + (ms->m1.t % MODULO_BASE)/2, ms->m2.t % MODULO_BASE);
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

int assert_float_is_equal_to(double estimate, double true_value, double abs_tolerance) {
    if ((estimate > (true_value - abs_tolerance)) && (estimate < (true_value + abs_tolerance))) {
        printf("\033[32mASSERTION PASSED:\033[0m Estimate lies within expected bounds from true value!\n");
        return 0; 
    } else {
        fprintf(stderr, "\033[31mASSERTION FAILED:\033[0m\n");
        fprintf(stderr, "ERROR: Estimate lies outside expected bounds from true value!\n");
        fprintf(stderr, "Expected bounds: %.5e...%.5e and received %.5e\n", true_value - abs_tolerance, true_value + abs_tolerance, estimate);
        return 1; 
    }

    UNREACHABLE("");
}


