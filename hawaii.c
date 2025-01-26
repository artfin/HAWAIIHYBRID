#include "hawaii.h"

#define HISTOGRAM_MAX_TPS 50

dipolePtr dipole = NULL;

MoleculeSystem init_ms(double mu, MonomerType t1, MonomerType t2, double *II1, double *II2, size_t seed) 
{
    INIT_WRANK;

    MoleculeSystem ms = {0};
    ms.mu = mu;

    ms.m1.t = t1;

    switch (t1) {
        case ATOM: break;
        case LINEAR_MOLECULE: {
          assert(II1[0] == II1[1]);
          memcpy(ms.m1.II, II1, 2*sizeof(double));
          break;
        }
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: {
          assert(II1[0] == II1[1]);
          memcpy(ms.m1.II, II1, 2*sizeof(double));
          break;
        }
        case ROTOR: {
          memcpy(ms.m1.II, II1, 3*sizeof(double));
          break;
        } 
        case ROTOR_REQUANTIZED_ROTATION: {
          TODO("init_ms");
        }
    } 
    
    ms.m1.qp   = malloc((t1%MODULO_BASE)   * sizeof(double));
    ms.m1.dVdq = malloc((t1%MODULO_BASE)/2 * sizeof(double));

    ms.m2.t = t2;
    switch (t2) {
        case ATOM: break;
        case LINEAR_MOLECULE: {
          assert(II2[0] == II2[1]);
          memcpy(ms.m2.II, II2, 2*sizeof(double));
          break;
        }
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: {
          assert(II2[0] == II2[1]);
          memcpy(ms.m2.II, II2, 2*sizeof(double));
          break;
        }
        case ROTOR: {
          memcpy(ms.m2.II, II2, 3*sizeof(double));
          break;
        } 
        case ROTOR_REQUANTIZED_ROTATION: {
          TODO("init_ms");
        }
    }

    ms.m2.qp   = malloc((t2%MODULO_BASE)   * sizeof(double));
    ms.m2.dVdq = malloc((t2%MODULO_BASE)/2 * sizeof(double));
    
    ms.QP_SIZE = (t1%MODULO_BASE) + (t2%MODULO_BASE) + 6;
    ms.Q_SIZE  = ms.QP_SIZE / 2;

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
        case LINEAR_MOLECULE:                      PRINT0("%.3e %.3e\n", ms.m1.II[0], ms.m1.II[1]); break;
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: PRINT0("%.3e %.3e\n", ms.m1.II[0], ms.m1.II[1]); break;
        case ROTOR:                                PRINT0("%.3e %.3e %.3e\n", ms.m1.II[0], ms.m1.II[1], ms.m1.II[2]); break;
        case ROTOR_REQUANTIZED_ROTATION:           PRINT0("%.3e %.3e %.3e\n", ms.m1.II[0], ms.m1.II[1], ms.m1.II[2]); break;

    }
    
    PRINT0("2nd monomer inertia tensor [%s]: ", monomer_type_name(t2));
    switch (t2) {
        case ATOM:                                 PRINT0("\n"); break;
        case LINEAR_MOLECULE:                      PRINT0("%.3e %.3e\n", ms.m2.II[0], ms.m2.II[1]); break;
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: PRINT0("%.3e %.3e\n", ms.m2.II[0], ms.m2.II[1]); break;
        case ROTOR:                                PRINT0("%.3e %.3e %.3e\n", ms.m2.II[0], ms.m2.II[1], ms.m2.II[2]); break;
        case ROTOR_REQUANTIZED_ROTATION:           PRINT0("%.3e %.3e %.3e\n", ms.m2.II[0], ms.m2.II[1], ms.m2.II[2]); break;
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

//void make_qp_odd(double *q, double *qp, size_t QP_SIZE) {
//    for (size_t k = 0; k < QP_SIZE; k += 2) {
//        qp[k + 1] = q[k / 2];
//    }
//} 

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

           deriv[IPHI]    = pPhi / m->II[0] / sin_theta / sin_theta;
           deriv[IPPHI]   = -m->dVdq[IPHI/2];
           deriv[ITHETA]  = pTheta / m->II[0]; 
           deriv[IPTHETA] = pPhi * pPhi * cos_theta / m->II[0] / sin_theta / sin_theta / sin_theta - m->dVdq[ITHETA/2]; 
          
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

           deriv[IPHI]    = pPhi / m->II[0] / sin_theta / sin_theta;
           deriv[IPPHI]   = -m->dVdq[IPHI/2];
           deriv[ITHETA]  = pTheta / m->II[0]; 
           deriv[IPTHETA] = pPhi * pPhi * cos_theta / m->II[0] / sin_theta / sin_theta / sin_theta - m->dVdq[ITHETA/2]; 
    
           break;                                                    
        }
        case ROTOR: {
          double pPhi   = m->qp[IPPHI];
          double theta  = m->qp[ITHETA];
          double pTheta = m->qp[IPTHETA];
          double psi    = m->qp[IPSI];
          double pPsi   = m->qp[IPPSI];

          double sin_theta = sin(theta);
          double cos_theta = cos(theta);

          double sin_psi = sin(psi); 
          double cos_psi = cos(psi);
          
          double t1 = (pPhi - pPsi * cos_theta) * sin_psi + pTheta * sin_theta * cos_psi;
          double t2 = (pPhi - pPsi * cos_theta) * cos_psi - pTheta * sin_theta * sin_psi; 
          deriv[IPHI] = t1 * sin_psi / (m->II[0] * sin_theta * sin_theta) + t2 * cos_psi / (m->II[1] * sin_theta * sin_theta);
          deriv[IPPHI] = -m->dVdq[IPHI/2]; 
   
          deriv[ITHETA] = t1 * cos_psi / (m->II[0] * sin_theta) - t2 * sin_psi / (m->II[1] * sin_theta);

          double t3 = pPhi * cos_theta - pPsi;
          double t4 = pPsi * cos_theta - pPhi; 
          deriv[IPTHETA] = -t3 * (t4 * (m->II[0] - m->II[1]) * cos_psi * cos_psi + \
                            pTheta * (m->II[0] - m->II[1]) * sin_theta * sin_psi * cos_psi + \
                            t4 * m->II[1]) / (m->II[0] * m->II[1] * sin_theta * sin_theta * sin_theta) - m->dVdq[ITHETA/2];

          deriv[IPSI] = -t1 * cos_theta * sin_psi / (m->II[0] * sin_theta * sin_theta) - \
                         t2 * cos_theta * cos_psi / (m->II[1] * sin_theta * sin_theta) + \
                         pPsi / m->II[2];
          deriv[IPPSI] = -t1*t2 / (m->II[0] * sin_theta * sin_theta) + t1*t2 / (m->II[1] * sin_theta * sin_theta) - m->dVdq[IPSI/2];
           
          break;                                                    
        }
        case ROTOR_REQUANTIZED_ROTATION: {
          TODO("rhsMonomer");
        }
    } 
}

void extract_dVdq_and_write_into_monomers(MoleculeSystem *ms) {
    memcpy(ms->m1.dVdq, ms->dVdq + 3,                            (ms->m1.t%MODULO_BASE)/2 * sizeof(double));
    memcpy(ms->m2.dVdq, ms->dVdq + 3 + (ms->m1.t%MODULO_BASE)/2, (ms->m2.t%MODULO_BASE)/2 * sizeof(double));
    
    // for (size_t i = 0; i < ms->m1.t%MODULO_BASE/2; ++i) {
    //   printf("m1.dVdq[%zu] = %.10e\n", i, ms->m1.dVdq[i]);
    // }

    // for (size_t i = 0; i < ms->m2.t%MODULO_BASE/2; ++i) {
    //   printf("m2.dVdq[%zu] = %.10e\n", i, ms->m2.dVdq[i]);
    // }
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
    
    double rhs_monomer1[ms->m1.t%MODULO_BASE];
    rhsMonomer(&ms->m1, rhs_monomer1);
    for (size_t i = 0; i < ms->m1.t%MODULO_BASE; ++i) {
        NV_Ith_S(ydot, i + 6) = rhs_monomer1[i];
    }

    double rhs_monomer2[ms->m2.t%MODULO_BASE];
    rhsMonomer(&ms->m2, rhs_monomer2);
    for (size_t i = 0; i < ms->m2.t%MODULO_BASE; ++i) {
        NV_Ith_S(ydot, i + 6 + (ms->m1.t%MODULO_BASE)) = rhs_monomer2[i];
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


gsl_matrix* compute_numerical_jac(void (*transform_angles)(double *qlab, double *qmol), double *qlab, size_t ninput_coordinates, size_t noutput_coordinates, size_t order)
{
    gsl_matrix *jac = gsl_matrix_alloc(ninput_coordinates, noutput_coordinates);
    for (size_t i = 0; i < ninput_coordinates; ++i) {
        for (size_t j = 0; j < noutput_coordinates; ++j) {
            gsl_matrix_set(jac, i, j, 0.0);
        }
    }

    double *result = (double*) malloc(noutput_coordinates * sizeof(double));
    double *temp = (double*) malloc(noutput_coordinates * sizeof(double));

    double step = 1e-3;

    for (size_t i = 0; i < ninput_coordinates; ++i) {
        memset(result, 0, noutput_coordinates*sizeof(double));

        switch (order) {
            case 2: {
                double c = qlab[i];
                qlab[i] = c + step;
                transform_angles(qlab, temp);
                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    result[j] += 1.0/2.0 * temp[j]; 
                }

                qlab[i] = c - step;
                transform_angles(qlab, temp);
                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    result[j] -= 1.0/2.0 * temp[j]; 
                }

                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    gsl_matrix_set(jac, i, j, result[j]/step);
                }
                
                qlab[i] = c;
                break;
            }
            case 4: {
                // 1/12, −2/3, 0, 2/3, −1/12	
                double c = qlab[i]; 
                qlab[i] = c - 2*step;
                transform_angles(qlab, temp);
                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    result[j] += 1.0/12.0 * temp[j]; 
                }
                
                qlab[i] = c - step;
                transform_angles(qlab, temp);
                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    result[j] -= 2.0/3.0 * temp[j]; 
                }
                
                qlab[i] = c + step;
                transform_angles(qlab, temp);
                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    result[j] += 2.0/3.0 * temp[j]; 
                }

                qlab[i] = c + 2*step;
                transform_angles(qlab, temp);
                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    result[j] -= 1.0/12.0 * temp[j]; 
                }
                
                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    gsl_matrix_set(jac, i, j, result[j]/step);
                }
                
                qlab[i] = c;
                break;
            }
            case 6: {
                // −1/60, 3/20, −3/4, 0, 3/4, −3/20, 1/60
                double c = qlab[i];

                qlab[i] = c - 3*step;
                transform_angles(qlab, temp);
                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    result[j] -= 1.0/60.0 * temp[j]; 
                }
                
                qlab[i] = c - 2*step;
                transform_angles(qlab, temp);
                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    result[j] += 3.0/20.0 * temp[j]; 
                }
                
                qlab[i] = c - step;
                transform_angles(qlab, temp);
                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    result[j] -= 3.0/4.0 * temp[j]; 
                }
                
                qlab[i] = c + step;
                transform_angles(qlab, temp);
                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    result[j] += 3.0/4.0 * temp[j]; 
                }
                
                qlab[i] = c + 2*step;
                transform_angles(qlab, temp);
                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    result[j] -= 3.0/20.0 * temp[j]; 
                }
                
                qlab[i] = c + 3*step;
                transform_angles(qlab, temp);
                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    result[j] += 1.0/60.0 * temp[j]; 
                }
                
                for (size_t j = 0; j < noutput_coordinates; ++j) {
                    gsl_matrix_set(jac, i, j, result[j]/step);
                }
                        
                qlab[i] = c;
                break;
            }
            default: {
                assert(0 && "ERROR: only order = 2, 4 and 6 are implemented"); 
            }
        }
    }

    free(result);
    free(temp);
    
    return jac;
}

Array compute_numerical_derivatives(double (*f)(double *q), double *q, size_t len, size_t order)
{
    Array derivatives = create_array(len);
    
    double step = 1e-6;
    if (order == 4) {
        step = 1e-4;
    } else if (order == 6) {
        step = 1e-3;
    }

    for (size_t i = 0; i < len; ++i) {
        switch (order) {
            case 2: {
              // -1/2, 0, 1/2
              double c = q[i]; 
              double r = 0.0;

              q[i] = c + step;
              r += 1.0/2.0 * f(q);

              q[i] = c - step;
              r -= 1.0/2.0 * f(q);

              q[i] = c;
              derivatives.data[i] = r/step; 
              break; 
            }
            case 4: {
                // 1/12, −2/3, 0, 2/3, −1/12	
                double c = q[i];
                double r = 0.0;

                q[i] = c - 2*step;
                r += 1.0/12.0 * f(q);

                q[i] = c - step;
                r -= 2.0/3.0 * f(q);

                q[i] = c + step;
                r += 2.0/3.0 * f(q);

                q[i] = c + 2*step;
                r -= 1.0/12.0 * f(q);

                q[i] = c;
                derivatives.data[i] = r/step; 
                break;
            }
            case 6: {
                // −1/60, 3/20, −3/4, 0, 3/4, −3/20, 1/60
                double c = q[i];
                double r = 0.0;

                q[i] = c - 3*step;
                r -= 1.0/60.0 * f(q);

                q[i] = c - 2*step;
                r += 3.0/20.0 * f(q);

                q[i] = c - step;
                r -= 3.0/4.0 * f(q);

                q[i] = c + step;
                r += 3.0/4.0 * f(q);

                q[i] = c + 2*step;
                r -= 3.0/20.0 * f(q);

                q[i] = c + 3*step;
                r += 1.0/60.0 * f(q);

                q[i] = c;
                derivatives.data[i] = r/step;
                break;
            }
            default: {
                assert(0 && "ERROR: only order = 2, 4 and 6 are implemented"); 
            }
        }
    }

    return derivatives;
}

Array compute_numerical_rhs(MoleculeSystem *ms, size_t order) 
{
    Array derivatives = create_array(ms->QP_SIZE);
    Array qp = create_array(ms->QP_SIZE);
   
    double step = 1e-6;
    if (order == 6) {
        step = 1e-3;
    }

    for (size_t i = 0; i < ms->QP_SIZE; ++i) {
        get_qp_from_ms(ms, &qp);

        if (i % 2 == 0) {
            switch (order) {
              case 2: {
                // −1/2, 0, 1/2
                double c = qp.data[i + 1];
                double r = 0.0;

                qp.data[i + 1] = c + step;
                put_qp_into_ms(ms, qp);
                r += 0.5 * Hamiltonian(ms);
                
                qp.data[i + 1] = c - step;
                put_qp_into_ms(ms, qp);
                r -= 0.5 * Hamiltonian(ms);
     
                qp.data[i + 1] = c; 
                put_qp_into_ms(ms, qp);

                derivatives.data[i] = r/step;
                break;
              }
              case 4: {
                // 1/12, −2/3, 0, 2/3, −1/12	
                double c = qp.data[i + 1];
                double r = 0.0;

                qp.data[i + 1] = c - 2*step;
                put_qp_into_ms(ms, qp);
                r += 1.0/12.0 * Hamiltonian(ms);
                
                qp.data[i + 1] = c - step;
                put_qp_into_ms(ms, qp);
                r -= 2.0/3.0 * Hamiltonian(ms);
                
                qp.data[i + 1] = c + step;
                put_qp_into_ms(ms, qp);
                r += 2.0/3.0 * Hamiltonian(ms);
                
                qp.data[i + 1] = c + 2*step;
                put_qp_into_ms(ms, qp);
                r -= 1.0/12.0 * Hamiltonian(ms);

                qp.data[i + 1] = c; 
                put_qp_into_ms(ms, qp);

                derivatives.data[i] = r/step; 
                break;
              }
              case 6: {
                // −1/60, 3/20, −3/4, 0, 3/4, −3/20, 1/60
                double c = qp.data[i + 1];
                double r = 0.0;

                qp.data[i + 1] = c - 3*step;
                put_qp_into_ms(ms, qp);
                r -= 1.0/60.0 * Hamiltonian(ms);
                
                qp.data[i + 1] = c - 2*step;
                put_qp_into_ms(ms, qp);
                r += 3.0/20.0 * Hamiltonian(ms);
                
                qp.data[i + 1] = c - step;
                put_qp_into_ms(ms, qp);
                r -= 3.0/4.0 * Hamiltonian(ms);
                
                qp.data[i + 1] = c + step;
                put_qp_into_ms(ms, qp);
                r += 3.0/4.0 * Hamiltonian(ms);
                
                qp.data[i + 1] = c + 2*step;
                put_qp_into_ms(ms, qp);
                r -= 3.0/20.0 * Hamiltonian(ms);
                
                qp.data[i + 1] = c + 3*step;
                put_qp_into_ms(ms, qp);
                r += 1.0/60.0 * Hamiltonian(ms);

                qp.data[i + 1] = c; 
                put_qp_into_ms(ms, qp);

                derivatives.data[i] = r/step; 
                break;
              }
              default: {
                assert(false);
              }
            }
        } else {
            switch (order) {
              case 2: {
                double c = qp.data[i - 1];
                double r = 0.0;

                qp.data[i - 1] = c + step;
                put_qp_into_ms(ms, qp);
                r += 0.5 * Hamiltonian(ms);

                qp.data[i - 1] = c - step;
                put_qp_into_ms(ms, qp);
                r -= 0.5 * Hamiltonian(ms);

                qp.data[i - 1] = c; 
                put_qp_into_ms(ms, qp);

                derivatives.data[i] = -r/step;
                break;
              }
              case 4: {
                double c = qp.data[i - 1];
                double r = 0.0;

                qp.data[i - 1] = c - 2*step;
                put_qp_into_ms(ms, qp);
                r += 1.0/12.0 * Hamiltonian(ms);
                
                qp.data[i - 1] = c - step;
                put_qp_into_ms(ms, qp);
                r -= 2.0/3.0 * Hamiltonian(ms);
                
                qp.data[i - 1] = c + step;
                put_qp_into_ms(ms, qp);
                r += 2.0/3.0 * Hamiltonian(ms);
                
                qp.data[i - 1] = c + 2*step;
                put_qp_into_ms(ms, qp);
                r -= 1.0/12.0 * Hamiltonian(ms);

                qp.data[i - 1] = c; 
                put_qp_into_ms(ms, qp);

                derivatives.data[i] = -r/step; 
                break;
              }
              case 6: {
                // −1/60, 3/20, −3/4, 0, 3/4, −3/20, 1/60
                double c = qp.data[i - 1];
                double r = 0.0;

                qp.data[i - 1] = c - 3*step;
                put_qp_into_ms(ms, qp);
                r -= 1.0/60.0 * Hamiltonian(ms);
                
                qp.data[i - 1] = c - 2*step;
                put_qp_into_ms(ms, qp);
                r += 3.0/20.0 * Hamiltonian(ms);
                
                qp.data[i - 1] = c - step;
                put_qp_into_ms(ms, qp);
                r -= 3.0/4.0 * Hamiltonian(ms);
                
                qp.data[i - 1] = c + step;
                put_qp_into_ms(ms, qp);
                r += 3.0/4.0 * Hamiltonian(ms);
                
                qp.data[i - 1] = c + 2*step;
                put_qp_into_ms(ms, qp);
                r -= 3.0/20.0 * Hamiltonian(ms);
                
                qp.data[i - 1] = c + 3*step;
                put_qp_into_ms(ms, qp);
                r += 1.0/60.0 * Hamiltonian(ms);

                qp.data[i - 1] = c; 
                put_qp_into_ms(ms, qp);

                derivatives.data[i] = -r/step; 
                break;
              }
              default: {
                assert(false);
              }
            }
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
          double phi1t    = ms->m1.qp[IPHI]; UNUSED(phi1t);
          double pphi1t   = ms->m1.qp[IPPHI];
          double theta1t  = ms->m1.qp[ITHETA];
          double ptheta1t = ms->m1.qp[IPTHETA];
          
          double sin_theta1t = sin(theta1t);
           
          KIN2 = ptheta1t * ptheta1t / (2.0 * ms->m1.II[0]) + pphi1t * pphi1t / (2.0 * ms->m1.II[1] * sin_theta1t * sin_theta1t);
          break;
        }
        case ROTOR_REQUANTIZED_ROTATION: 
        case ROTOR: {
          double phi1t    = ms->m1.qp[IPHI]; UNUSED(phi1t);
          double pphi1t   = ms->m1.qp[IPPHI];
          double theta    = ms->m1.qp[ITHETA];
          double ptheta1t = ms->m1.qp[IPTHETA];
          double psi      = ms->m1.qp[IPSI];
          double ppsi1t   = ms->m1.qp[IPPSI];

          double sin_theta = sin(theta);
          double cos_theta = cos(theta);

          double sin_psi = sin(psi); 
          double cos_psi = cos(psi);

          KIN2 = ((pphi1t - ppsi1t * cos_theta) * sin_psi + ptheta1t * sin_theta * cos_psi) * \
                 ((pphi1t - ppsi1t * cos_theta) * sin_psi + ptheta1t * sin_theta * cos_psi) / (2.0 * ms->m1.II[0] * sin_theta * sin_theta) + \
                 ((pphi1t - ppsi1t * cos_theta) * cos_psi - ptheta1t * sin_theta * sin_psi) * \
                 ((pphi1t - ppsi1t * cos_theta) * cos_psi - ptheta1t * sin_theta * sin_psi) / (2.0 * ms->m1.II[1] * sin_theta * sin_theta) + \
                 ppsi1t * ppsi1t / (2.0 * ms->m1.II[2]);
          break;
        }
    }
    
    switch (ms->m2.t) {
        case ATOM: break;
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: 
        case LINEAR_MOLECULE: { 
          double phi2t    = ms->m2.qp[IPHI]; UNUSED(phi2t);
          double pphi2t   = ms->m2.qp[IPPHI];
          double theta2t  = ms->m2.qp[ITHETA];
          double ptheta2t = ms->m2.qp[IPTHETA];
          
          double sin_theta2t = sin(theta2t);
    
          KIN3 = ptheta2t * ptheta2t / (2.0 * ms->m2.II[0]) + pphi2t * pphi2t / (2.0 * ms->m2.II[1] * sin_theta2t * sin_theta2t);
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
    extract_q(ms->m2.qp,             ms->intermediate_q+6/2+(ms->m1.t%MODULO_BASE)/2, ms->m2.t%MODULO_BASE);
}

void invert_momenta(MoleculeSystem *ms) {
    for (size_t i = 1; i < 6; i += 2) {
        ms->intermolecular_qp[i] = -ms->intermolecular_qp[i];
    }

    for (size_t i = 1; i < (ms->m1.t%MODULO_BASE); i += 2) {
        ms->m1.qp[i] = -ms->m1.qp[i];
    }
    
    for (size_t i = 1; i < (ms->m2.t%MODULO_BASE); i += 2) {
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
    bool success = false;

    double r;

    do {
        double U = mt_drand();
        double V = mt_drand();

        double logU = log(U);
        if (isnan(logU) || (isinf(logU) != 0)) {
            // printf("WARNING: generate U failed\n"); 
            continue;    
        }

        r = sigma * sqrt(-2 * logU) * cos(2.0 * M_PI * V);
        success = true;
    } while(!success);

    return r;
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
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: {
          TODO("q_generator");
        }
        case ROTOR: {
          ms->m1.qp[IPHI]   = mt_drand() * 2.0 * M_PI;
          ms->m1.qp[ITHETA] = acos(2.0*mt_drand() - 1.0);
          ms->m1.qp[IPSI]   = mt_drand() * 2.0 * M_PI;
          break;
        }                        
        case ROTOR_REQUANTIZED_ROTATION: { 
          TODO("q_generator");
        }
    } 
    
    switch (ms->m2.t) {
        case ATOM: break;
        case LINEAR_MOLECULE: {
          ms->m2.qp[IPHI]   = mt_drand() * 2.0 * M_PI;
          ms->m2.qp[ITHETA] = acos(2.0*mt_drand() - 1.0);
          break;
        }
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: {
          TODO("q_generator");
        }
        case ROTOR: {
          ms->m2.qp[IPHI]   = mt_drand() * 2.0 * M_PI;
          ms->m2.qp[ITHETA] = acos(2.0*mt_drand() - 1.0);
          ms->m2.qp[IPSI]   = mt_drand() * 2.0 * M_PI;
          break;
        }                        
        case ROTOR_REQUANTIZED_ROTATION: { 
          TODO("q_generator");
        }
    }
}

static void p_generator_linear_molecule(Monomer *m, double Temperature)
{
    double sqrt_IIKT = sqrt(m->II[0] / HkT * Temperature);
    double sin_theta = sin(m->qp[ITHETA]); 

    double x0 = generate_normal(1.0);
    double x1 = generate_normal(1.0);
   
    assert(m->II[0] == m->II[1]);
    m->qp[IPPHI]   = sqrt_IIKT * sin_theta * x1;
    m->qp[IPTHETA] = sqrt_IIKT * x0;
}

static void p_generator_rotor(Monomer *m, double Temperature)
{
    double theta = m->qp[ITHETA]; 
    double psi = m->qp[IPSI]; 
    
    double sqrt_IIx = sqrt(m->II[0] / HkT * Temperature);
    double sqrt_IIy = sqrt(m->II[1] / HkT * Temperature);
    double sqrt_IIz = sqrt(m->II[2] / HkT * Temperature);

    double sin_theta, cos_theta;
    double sin_psi, cos_psi;

    sincos(theta, &sin_theta, &cos_theta);
    sincos(psi, &sin_psi, &cos_psi);
    
    double x0 = generate_normal(1.0);
    double x1 = generate_normal(1.0);
    double x2 = generate_normal(1.0);

    m->qp[IPPHI]   = sqrt_IIx * sin_theta * sin_psi * x0 + sqrt_IIy * sin_theta * cos_psi * x1 + sqrt_IIz * cos_theta * x2;
    m->qp[IPTHETA] = sqrt_IIx * cos_psi * x0 - sqrt_IIy * sin_psi * x1;
    m->qp[IPPSI]   = sqrt_IIz * x2;
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
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: {
          TODO("p_generator");
        }
        case ROTOR: {
          p_generator_rotor(&ms->m1, Temperature);
          break;
        }
        case ROTOR_REQUANTIZED_ROTATION: { 
          TODO("p_generator");
        }
    }

    switch (ms->m2.t) {
        case ATOM: break;
        case LINEAR_MOLECULE: {
          p_generator_linear_molecule(&ms->m2, Temperature);
          break;
        }
        case LINEAR_MOLECULE_REQUANTIZED_ROTATION: {
          TODO("p_generator");
        }
        case ROTOR: {
          p_generator_rotor(&ms->m2, Temperature);
          break;
        }
        case ROTOR_REQUANTIZED_ROTATION: { 
          TODO("p_generator");
        }
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

void compute_dHdp(MoleculeSystem *ms, gsl_matrix* dHdp) 
{
    double R      = ms->intermolecular_qp[IR];
    double pR     = ms->intermolecular_qp[IPR];
    double pPhi   = ms->intermolecular_qp[IPPHI];
    double Theta  = ms->intermolecular_qp[ITHETA];
    double pTheta = ms->intermolecular_qp[IPTHETA];
    
    double sinTheta = sin(Theta);

    // strange choice to use /2... 
    gsl_matrix_set(dHdp, 0, IPHI/2,   pPhi / (ms->mu * R * R * sinTheta * sinTheta));
    gsl_matrix_set(dHdp, 0, ITHETA/2, pTheta / (ms->mu * R * R));
    gsl_matrix_set(dHdp, 0, IR/2,     pR / ms->mu);
    
    double rhs_monomer1[ms->m1.t%MODULO_BASE];
    rhsMonomer(&ms->m1, rhs_monomer1);
    for (size_t i = 0; i < ms->m1.t%MODULO_BASE; i += 2) {
        gsl_matrix_set(dHdp, 0, i/2 + 6/2, rhs_monomer1[i]);
    }

    double rhs_monomer2[ms->m2.t%MODULO_BASE];
    rhsMonomer(&ms->m2, rhs_monomer2);
    for (size_t i = 0; i < ms->m2.t%MODULO_BASE; i += 2) {
        gsl_matrix_set(dHdp, 0, i/2 + 6/2 + ms->m1.t%MODULO_BASE/2, rhs_monomer2[i]);
    }
}

void calculate_M2(MoleculeSystem *ms, CalcParams *params, double Temperature, double *m, double *q)
// Running mean/variance formulas taken from GSL 1.15
// https://github.com/ampl/gsl/blob/master/monte/plain.c 
{
    assert(params->initialM2_npoints > 0);
    assert(fabs(params->pesmin) > 1e-15);

    double h = 1.0e-3;

    gsl_matrix *D           = gsl_matrix_alloc(ms->Q_SIZE, 3);
    gsl_matrix *dHdp        = gsl_matrix_alloc(1, ms->Q_SIZE);
    gsl_matrix *dip_lab_dot = gsl_matrix_alloc(1, 3);
    
    size_t counter = 0;
    size_t desired_dist = 0;
    size_t integral_counter = 0;

    *m = 0.0;
    *q = 0.0;

    size_t print_every_nth_iteration = 1;
    switch (params->ps) {
        case FREE_AND_METASTABLE: print_every_nth_iteration = 100000; break;
        case BOUND:               print_every_nth_iteration = 1000; break;
    } 

    while (integral_counter < params->initialM2_npoints) {
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
    
            double dp[3];
            double dm[3];

            for (size_t i = 0; i < ms->Q_SIZE; ++i) {
                double tmp = ms->intermediate_q[i];
                
                ms->intermediate_q[i] = tmp + h;
                (*dipole)(ms->intermediate_q, dp);
                
                ms->intermediate_q[i] = tmp - h;
                (*dipole)(ms->intermediate_q, dm);

                gsl_matrix_set(D, i, 0, (dp[0] - dm[0])/(2.0*h)); 
                gsl_matrix_set(D, i, 1, (dp[1] - dm[1])/(2.0*h)); 
                gsl_matrix_set(D, i, 2, (dp[2] - dm[2])/(2.0*h)); 

                ms->intermediate_q[i] = tmp;
            } 

            compute_dHdp(ms, dHdp); 

            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, dHdp, D, 0.0, dip_lab_dot);
            
            double fval = gsl_matrix_get(dip_lab_dot, 0, 0)*gsl_matrix_get(dip_lab_dot, 0, 0) + \
                          gsl_matrix_get(dip_lab_dot, 0, 1)*gsl_matrix_get(dip_lab_dot, 0, 1) + \
                          gsl_matrix_get(dip_lab_dot, 0, 2)*gsl_matrix_get(dip_lab_dot, 0, 2); 
            double diff = fval - *m;
            *m += diff / (integral_counter + 1.0);
            *q += diff * diff * (integral_counter / (integral_counter + 1.0));
            integral_counter++;

            if (integral_counter % print_every_nth_iteration == 0) {
                double M2_est = *m * SecondCoeff * params->partial_partition_function_ratio;
                double M2std_est = sqrt(*q / integral_counter / (integral_counter - 1)) * SecondCoeff * params->partial_partition_function_ratio;
                printf("[calculate_M2] %zu/%zu points: \t M2 = %.5e +/- %.5e\n", integral_counter, params->initialM2_npoints, M2_est, M2std_est);
            }
        }
    } 

    *m = *m * SecondCoeff * params->partial_partition_function_ratio;
    *q = sqrt(*q / integral_counter / (integral_counter - 1)) * SecondCoeff * params->partial_partition_function_ratio;
}

#ifdef USE_MPI
void mpi_calculate_M0(MoleculeSystem *ms, CalcParams *params, double Temperature, double *m, double *q)
{
    assert(params->initialM0_npoints > 0);
    assert(fabs(params->pesmin) > 1e-15);
    assert(params->partial_partition_function_ratio > 0);

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

void mpi_calculate_M2(MoleculeSystem *ms, CalcParams *params, double Temperature, double *m, double *q)
{
    assert(params->initialM2_npoints > 0);
    assert(fabs(params->pesmin) > 1e-15);
    assert(params->partial_partition_function_ratio > 0);

    INIT_WRANK;
    INIT_WSIZE;
    
    double h = 1.0e-3;

    gsl_matrix *D           = gsl_matrix_alloc(ms->Q_SIZE, 3);
    gsl_matrix *dHdp        = gsl_matrix_alloc(1, ms->Q_SIZE);
    gsl_matrix *dip_lab_dot = gsl_matrix_alloc(1, 3);

    size_t counter = 0;
    size_t desired_dist = 0;
    size_t integral_counter = 0;

    size_t print_every_nth_iteration = 1;
    switch (params->ps) {
        case FREE_AND_METASTABLE: print_every_nth_iteration = 1000000; break;
        case BOUND:               print_every_nth_iteration = 1000;    break;
    } 
    
    size_t local_npoints = params->initialM2_npoints / _wsize;
    
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
            
            double dp[3];
            double dm[3];

            for (size_t i = 0; i < ms->Q_SIZE; ++i) {
                double tmp = ms->intermediate_q[i];
                
                ms->intermediate_q[i] = tmp + h;
                (*dipole)(ms->intermediate_q, dp);
                
                ms->intermediate_q[i] = tmp - h;
                (*dipole)(ms->intermediate_q, dm);

                gsl_matrix_set(D, i, 0, (dp[0] - dm[0])/(2.0*h)); 
                gsl_matrix_set(D, i, 1, (dp[1] - dm[1])/(2.0*h)); 
                gsl_matrix_set(D, i, 2, (dp[2] - dm[2])/(2.0*h)); 

                ms->intermediate_q[i] = tmp;
            } 

            compute_dHdp(ms, dHdp); 

            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, dHdp, D, 0.0, dip_lab_dot);
            
            double fval = gsl_matrix_get(dip_lab_dot, 0, 0)*gsl_matrix_get(dip_lab_dot, 0, 0) + \
                          gsl_matrix_get(dip_lab_dot, 0, 1)*gsl_matrix_get(dip_lab_dot, 0, 1) + \
                          gsl_matrix_get(dip_lab_dot, 0, 2)*gsl_matrix_get(dip_lab_dot, 0, 2);

            if (isnan(fval) || (isinf(fval) != 0)) {
                PRINT0("Energy = %.10e\n", energy);
                PRINT0("Caught NaN value of integrand at:\n");
                
                Array qp = create_array(ms->QP_SIZE);
                get_qp_from_ms(ms, &qp);

                for (size_t i = 0; i < ms->QP_SIZE; ++i) {
                    PRINT0("qp[%zu] = %.10e\n", i, qp.data[i]);
                }

                free_array(&qp);
                    
                continue;
            }

            double diff = fval - ml;
            ml += diff / (integral_counter + 1.0);
            ql += diff * diff * (integral_counter / (integral_counter + 1.0));
            integral_counter++;

            if (integral_counter % print_every_nth_iteration == 0) {
                double M2_est    = ml * SecondCoeff * params->partial_partition_function_ratio;
                double M2std_est = sqrt(ql / integral_counter / (integral_counter - 1)) * SecondCoeff * params->partial_partition_function_ratio;
                PRINT0("[mpi_calculate_M2] %zu/%zu points: \t M2 = %.5e +/- %.5e\n", _wsize*integral_counter, params->initialM2_npoints, M2_est, M2std_est);
            }
        }
    } 

    MPI_Allreduce(MPI_IN_PLACE, &ml, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &ql, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
     
    *m = (ml/_wsize) * SecondCoeff * params->partial_partition_function_ratio;
    *q = sqrt((ql/_wsize) / integral_counter / (integral_counter - 1)) * SecondCoeff * params->partial_partition_function_ratio;
    
    gsl_matrix_free(D);
    gsl_matrix_free(dHdp);
    gsl_matrix_free(dip_lab_dot);
}
#endif // USE_MPI


#include "trajectory.h"

void track_turning_points(Tracker *tr, double R)
{
    tr->before2 = tr->before;
    tr->before  = tr->current;
    tr->current = R; 

    if (tr->ready) {
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
    if (tr->called > 1) tr->ready = true; 
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
   
    // return value: 
    // 0 -- means trajectories propagated successfully
    // >0 -- an error was encountered during trajectory propagation 
    int status = 0;
    
    double t = 0.0;
    double tout = params->sampling_time;
   
    Tracker tr = {
      .before2 = qp.data[IR],
      .before  = qp.data[IR],
      .current = qp.data[IR],
      .ready   = false,
    };

    double prev_value, curr_value;

    /*
     * We start step_counter from 1 so that correlation value after the first integration step
     * will go into correlation_forw[1] 
     */
    for (size_t step_counter = 1; step_counter < params->MaxTrajectoryLength; ++step_counter, tout += params->sampling_time)
    {
        status = make_step(traj, tout, &t);
        if (status) {
            printf("CVODE ERROR: status = %d\n", status);
            break;
        }

        put_qp_into_ms(ms, (Array){.data = N_VGetArrayPointer(traj->y), .n = ms->QP_SIZE});
        extract_q_and_write_into_ms(ms);
        (*dipole)(ms->intermediate_q, dipt);
        
        if (isnan(dipt[0]) || isnan(dipt[1]) || isnan(dipt[2])) {
            printf("ERROR: one of the components of the dipole is corrupted!\n");
            printf("The initial phase-point for broken trajectory in the forward direction is:\n");
            for (size_t i = 0; i < ms->QP_SIZE; ++i) {
                printf("%.10e ", qp.data[i]);
            }
            printf("\n");
            return 1;         
        }
        
        prev_value = curr_value;
        curr_value = dip0[0]*dipt[0] + dip0[1]*dipt[1] + dip0[2]*dipt[2];

        if (step_counter > 1) {
            double ratio = fabs(curr_value / prev_value);
            if (ratio > 10000) {
                printf("ratio = %.10f, prev = %.10e, curr = %.10e\n", ratio, prev_value, curr_value); 
            }
            if (ratio > 1e10) {
                printf("ERROR: unexpectedly large jump in dipole value!\n"); 
                printf("The initial phase-point for broken trajectory in the forward direction is:\n");
                for (size_t i = 0; i < ms->QP_SIZE; ++i) {
                    printf("%.10e ", qp.data[i]);
                }
                printf("\n");
                return 1;         
            } 
        }

        correlation_forw[step_counter] = curr_value; 

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
    tr.ready   = false;

    for (size_t step_counter = 1; step_counter < params->MaxTrajectoryLength; ++step_counter, tout += params->sampling_time)
    {
        status = make_step(traj, tout, &t);
        if (status) {
            printf("CVODE ERROR: status = %d\n", status);
            break;
        }

        put_qp_into_ms(ms, (Array){.data = N_VGetArrayPointer(traj->y), .n = ms->QP_SIZE});
        extract_q_and_write_into_ms(ms);
        (*dipole)(ms->intermediate_q, dipt);
        
        if (isnan(dipt[0]) || isnan(dipt[1]) || isnan(dipt[2])) {
            printf("ERROR: one of the components of the dipole is corrupted!\n");
            printf("The initial phase-point for broken trajectory in the backward direction is:\n");
            for (size_t i = 0; i < ms->QP_SIZE; ++i) {
                printf("%.10e ", qp.data[i]);
            }
            printf("\n");
            return 1;         
        }
        
        prev_value = curr_value;
        curr_value = dip0[0]*dipt[0] + dip0[1]*dipt[1] + dip0[2]*dipt[2];

        if (step_counter > 1) {
            double ratio = fabs(curr_value / prev_value); 
            if (ratio > 10000) {
                printf("ratio = %.10f, prev = %.10e, curr = %.10e\n", ratio, prev_value, curr_value); 
            }
            if (ratio > 1e10) {
                printf("ERROR: unexpectedly large jump in dipole value!\n"); 
                printf("The initial phase-point for broken trajectory in the backward direction is:\n");
                for (size_t i = 0; i < ms->QP_SIZE; ++i) {
                    printf("%.10e ", qp.data[i]);
                }
                printf("\n");
                return 1;         
            } 
        }


        correlation_back[step_counter] = curr_value; 
        
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
  
    FILE *fp = NULL; 
    if (_wrank == 0) {
        fp = fopen(params->cf_filename, "w");
        if (fp == NULL) { 
            printf("ERROR: Could not open '%s' for writing! Exiting...\n", params->cf_filename);
            exit(1);
        }
    } 

    // TODO: handle the distribution of trajectories between processes so that in total we 
    //       calculate 'total_trajectories/niterations' each iteration. Basically, handle 
    //       the case when total_trajectories is not divisible by niterations. 
    size_t local_ntrajectories = params->total_trajectories / params->niterations / _wsize;
   
    double *crln       = malloc(params->MaxTrajectoryLength * sizeof(double));
    double *local_crln = malloc(params->MaxTrajectoryLength * sizeof(double));
    
    CFnc total_crln = {
        .t           = linspace(0.0, params->sampling_time*(params->MaxTrajectoryLength-1), params->MaxTrajectoryLength),
        .data        = malloc(params->MaxTrajectoryLength * sizeof(double)),
        .len         = params->MaxTrajectoryLength,
        .ntraj       = 0,
        .Temperature = Temperature,
    }; 
    memset(total_crln.data, 0, params->MaxTrajectoryLength * sizeof(double)); 
    
    CFnc total_crln_iter = {
        .t           = linspace(0.0, params->sampling_time*(params->MaxTrajectoryLength-1), params->MaxTrajectoryLength),
        .data        = malloc(params->MaxTrajectoryLength * sizeof(double)),
        .len         = params->MaxTrajectoryLength,
        .ntraj       = 0,
        .Temperature = Temperature,
    }; 
    memset(total_crln_iter.data, 0, params->MaxTrajectoryLength * sizeof(double)); 
    
    Trajectory traj = init_trajectory(ms, params->cvode_tolerance);

    // TODO: this histogram is filled but not extended yet at the appropriate moment
    // TODO: we should probably save it to a file at the end of the iteration if 
    //       some variable is set in the 'params' structure.
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
    PRINT0("Running preliminary calculations of M0 & M2 using rejection sampler to generate phase-points from Boltzmann distribution\n");
    PRINT0("The estimate for M0 will be based on %zu points\n", params->initialM0_npoints); 
    PRINT0("The estimate for M2 will be based on %zu points\n\n", params->initialM2_npoints); 

    double prelim_M0, prelim_M0std;
    mpi_calculate_M0(ms, params, Temperature, &prelim_M0, &prelim_M0std);
    PRINT0("M0 = %.10e +/- %.10e [%.10e ... %.10e]\n", prelim_M0, prelim_M0std, prelim_M0 - prelim_M0std, prelim_M0 + prelim_M0std);
    PRINT0("Error: %.3f%%\n", prelim_M0std/prelim_M0 * 100.0);
    
    double prelim_M2, prelim_M2std;
    mpi_calculate_M2(ms, params, Temperature, &prelim_M2, &prelim_M2std); 
    PRINT0("M2 = %.10e +/- %.10e [%.10e ... %.10e]\n", prelim_M2, prelim_M2std, prelim_M2-prelim_M2std, prelim_M2+prelim_M2std);
    PRINT0("Error: %.3f%%\n", prelim_M2std/prelim_M2 * 100.0);
    
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
            }
        }

        MPI_Allreduce(local_crln, total_crln_iter.data, params->MaxTrajectoryLength, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
       
        for (size_t i = 0; i < params->MaxTrajectoryLength; ++i) {
            total_crln.data[i] += total_crln_iter.data[i];
        }
        total_crln.ntraj += local_ntrajectories * _wsize;

        memset(local_crln,           0, params->MaxTrajectoryLength * sizeof(double));
        memset(total_crln_iter.data, 0, params->MaxTrajectoryLength * sizeof(double));

        PRINT0("ITERATION %zu/%zu: accumulated %zu trajectories. Saving the temporary result to '%s'\n", iter+1, params->niterations, total_crln.ntraj, params->cf_filename);
        double M0_crln_est = total_crln.data[0] / total_crln.ntraj * ZeroCoeff;
        PRINT0("M0 ESTIMATE FROM CF: %.5e, PRELIMINARY M0 ESTIMATE: %.5e, diff: %.3f%%\n", M0_crln_est, prelim_M0, (M0_crln_est - prelim_M0)/prelim_M0*100.0);

        // double M2_crln_est = SecondCoeff * 2.0/params->sampling_time/params->sampling_time*(total_crln.data[0] - total_crln.data[1]);
        double M2_crln_est = -SecondCoeff * (35.0*total_crln.data[0] - 104.0*total_crln.data[1] + 114.0*total_crln.data[2] - 56.0*total_crln.data[3] + 11.0*total_crln.data[4])/12.0/params->sampling_time/params->sampling_time/1000;

        PRINT0("M2 ESTIMATE FROM CF: %.5e, PRELIMINARY M2 ESTIMATE: %.5e, diff: %.3f%%\n\n", M2_crln_est, prelim_M2, (M2_crln_est - prelim_M2)/prelim_M2*100.0);


        if (_wrank == 0) { 
            save_correlation_function(fp, total_crln, params);
        }
    }
 
    if (tps_hist != NULL) {
        gsl_histogram_free(tps_hist);
    }
    
    free(crln);
    free(local_crln);
    free_cfnc(total_crln_iter);

    for (size_t i = 0; i < params->MaxTrajectoryLength; ++i) {
        total_crln.data[i] /= total_crln.ntraj;
    }

    return total_crln; 
}
#endif // USE_MPI

void gsl_fft_square(double *farr, size_t N) {
    farr[0] = farr[0] * farr[0];
    farr[N/2] = farr[N/2] * farr[N/2];

    for (size_t i = 1; i < N/2; ++i) {
        farr[i] = farr[i] * farr[i] + farr[N-i] * farr[N-i];
    }
}

#ifdef USE_MPI
SFnc calculate_spectral_function_using_prmu_representation_and_save(MoleculeSystem *ms, CalcParams *params, double Temperature) 
{
    // TODO: should we do any M0/M2 estimates at the beginning and then show the convergence throughtout the iterations? 
    assert(dipole != NULL);

    assert(params->MaxTrajectoryLength > 0);
    assert(params->sampling_time > 0);
    assert(params->cvode_tolerance > 0);
    assert(params->total_trajectories > 0);
    assert(params->niterations >= 1);
    assert(params->ApproximateFrequencyMax > 0);
    assert(params->sf_filename != NULL);
    
    INIT_WRANK;
    INIT_WSIZE;
    
    FILE *fp = NULL; 
    if (_wrank == 0) {
        fp = fopen(params->sf_filename, "w");
        if (fp == NULL) { 
            printf("ERROR: Could not open '%s' for writing! Exiting...\n", params->sf_filename);
            exit(1);
        }
    } 

    double frequency_step = 1.0 / (params->sampling_time * ATU) / LightSpeed_cm / params->MaxTrajectoryLength; // cm^-1
    double theoretical_frequency_max = 1.0 / (LightSpeed * 100.0 * params->sampling_time * ATU) / 2.0; // cm^-1

    if (params->ApproximateFrequencyMax > theoretical_frequency_max) {
        PRINT0("ERROR: Requested maximum frequency (%.3e) should be less than 1/2 of the maximum signal frequency, which is (%.3e)\n",
                params->ApproximateFrequencyMax, theoretical_frequency_max);
        exit(1);
    }

    size_t frequency_array_length = (int) (params->ApproximateFrequencyMax / frequency_step) + 1; 
    double max_frequency = frequency_step * (frequency_array_length - 1);

    PRINT0("\n\n");
    PRINT0("------------------------------------------------------------------------\n");
    PRINT0("Calculating spectral function using pr/mu representation at T = %.2f using following parameters:\n", Temperature);
    PRINT0("    trajectories to be calculated (total_trajectories):                  %zu\n",  params->total_trajectories);
    PRINT0("    # of iterations that the calculation is divided into (niterations):  %zu\n",  params->niterations);
    PRINT0("    maximum length of trajectory (MaxTrajectoryLength):                  %zu\n",  params->MaxTrajectoryLength);
    PRINT0("    sampling time of dipole on trajectory (sampling_time):               %.2f\n", params->sampling_time);
    PRINT0("    initial intermolecular distance (R0):                                %.2f\n", params->R0);
    PRINT0("    CVode tolerance:                                                     %.3e\n", params->cvode_tolerance);
    PRINT0("    Approximate maximum frequency:                                       %.3e\n", params->ApproximateFrequencyMax);
    PRINT0("------------------------------------------------------------------------\n");
    PRINT0("\n");

    PRINT0("Theoretical maximum frequency with provided parameters: %.3e\n", theoretical_frequency_max);
    PRINT0("The frequency step with provided parameters:            %.3e\n", frequency_step);
    PRINT0("The spectral function will be calculated in the range [0 .. %.3e] for %zu points\n", max_frequency, frequency_array_length);

    PRINT0("NOTE: Connes apodization will be applied to the time dependence of dipole before applying Fourier transform\n"); 

    // convert spectral funciton [atomic units -> J * m^6 * s]
    double SF_COEFF = Hartree * pow(ALU, 6) * ATU * params->R0 * params->R0 * params->sampling_time * params->sampling_time;

    size_t local_ntrajectories = params->total_trajectories / params->niterations / _wsize;

    double dipt[3]; 

    double *dipx = (double*) malloc(params->MaxTrajectoryLength * sizeof(double));
    assert(dipx != NULL && "ASSERT: not enough memory!\n");
    memset(dipx, 0, params->MaxTrajectoryLength * sizeof(double)); 

    double *dipy = (double*) malloc(params->MaxTrajectoryLength * sizeof(double));
    assert(dipy != NULL && "ASSERT: not enough memory!\n");
    memset(dipy, 0, params->MaxTrajectoryLength * sizeof(double)); 

    double *dipz = (double*) malloc(params->MaxTrajectoryLength * sizeof(double));
    assert(dipz != NULL && "ASSERT: not enough memory!\n"); 
    memset(dipz, 0, params->MaxTrajectoryLength * sizeof(double)); 

    SFnc sf_iter = {
        .nu       = NULL,
        .data     = (double*) malloc(frequency_array_length * sizeof(double)),
        .len      = frequency_array_length,
        .capacity = frequency_array_length,
        .ntraj    = 0,
    };
    assert(sf_iter.data != NULL && "ASSERT: not enough memory!\n"); 

    SFnc sf_total = {
        .nu       = linspace(0.0, max_frequency, frequency_array_length),
        .data     = (double*) malloc(frequency_array_length * sizeof(double)),
        .len      = frequency_array_length,
        .capacity = frequency_array_length,
        .ntraj    = 0,
    };

    assert(sf_total.nu != NULL && "ASSERT: not enough memory!\n"); 
    assert(sf_total.data != NULL && "ASSERT: not enough memory!\n"); 
    memset(sf_total.data, 0, frequency_array_length * sizeof(double)); 
    

    Trajectory traj = init_trajectory(ms, params->cvode_tolerance);
    
    int status = 0; 
    double t, tout;

    Array qp0 = create_array(ms->QP_SIZE);
   
    PRINT0("\n\n"); 
    PRINT0("Running preliminary calculation of M0 to test the sampler and dipole function...\n");
    PRINT0("The estimate will be based on %zu points\n\n", params->initialM0_npoints); 
    
    double prelim_M0, prelim_M0std;
    mpi_calculate_M0(ms, params, Temperature, &prelim_M0, &prelim_M0std);
    PRINT0("M0 = %.10e +/- %.10e [%.10e ... %.10e]\n", prelim_M0, prelim_M0std, prelim_M0 - prelim_M0std, prelim_M0 + prelim_M0std);
    PRINT0("Error: %.3f%%\n\n", prelim_M0std/prelim_M0 * 100.0);
    
    double prelim_M2, prelim_M2std; 
    mpi_calculate_M2(ms, params, Temperature, &prelim_M2, &prelim_M2std);
    PRINT0("M2 = %.10e +/- %.10e [%.10e ... %.10e]\n", prelim_M2, prelim_M2std, prelim_M2 - prelim_M2std, prelim_M2 + prelim_M2std);
    PRINT0("Error: %.3f%%\n\n", prelim_M2std/prelim_M2 * 100.0);

    for (size_t iter = 0; iter < params->niterations; ++iter) {

        memset(sf_iter.data, 0, frequency_array_length * sizeof(double));

        for (size_t traj_counter = 0; traj_counter < local_ntrajectories; ) 
        {
            q_generator(ms, params);
            ms->intermolecular_qp[IR] = params->R0;
            p_generator(ms, Temperature);

            if (ms->intermolecular_qp[IPR] > 0.0) {
                ms->intermolecular_qp[IPR] = -ms->intermolecular_qp[IPR]; 
            }

            double pr_mu = -ms->intermolecular_qp[IPR] / ms->mu;
            
            get_qp_from_ms(ms, &qp0);
            set_initial_condition(&traj, qp0);

            t    = 0.0;
            tout = params->sampling_time;

            for (size_t step_counter = 0; step_counter < params->MaxTrajectoryLength; ++step_counter, tout += params->sampling_time)
            {
                status = make_step(&traj, tout, &t);
                if (status) {
                    printf("CVODE ERROR: status =  %d\n", status);
                    break;
                }

                put_qp_into_ms(ms, (Array){
                    .data = N_VGetArrayPointer(traj.y),
                    .n    = ms->QP_SIZE,
                });
                extract_q_and_write_into_ms(ms);
                (*dipole)(ms->intermediate_q, dipt);

                dipx[step_counter] = dipt[0];
                dipy[step_counter] = dipt[1];
                dipz[step_counter] = dipt[2];

                // printf("t = %.2f, R = %.3e\n", t, ms->intermolecular_qp[IR]);
                // printf("t = %.2f, R = %.3e, dipx = %.3e, dipy = %.3e, dipz = %.3e\n", t, ms->intermolecular_qp[IR], dipx[step_counter], dipy[step_counter], dipz[step_counter]); 

                if (ms->intermolecular_qp[IR] > params->R0) {
                    break;
                }
            }

            if (status) {
                printf("Caught CVode error. Resampling new conditions...\n");
                continue;
            }
           
            // if we hit this 'if', it means that the length of the trajectory turned out to be longer
            // than MaxTrajectoryLength. Not sure if it's a good idea to add this trajectory to the 
            // result, but it probably does not matter, because for reasonably large MaxTrajectoryLength        
            // there will be only a few of those long-living metastable trajectories.  
            if (ms->intermolecular_qp[IR] < params->R0) {
                printf("INFO: the trajectory is too long. Fouier transform of the dipole from this trajectory will not be added to the cumulative result. Continuing...\n"); 
                continue;
            }

            // Here we chose to apply apodization procedure to dipole time-dependence
            // TODO: should we make it customizable? 
            connes_apodization((Array){.data = dipx, .n = params->MaxTrajectoryLength}, params->sampling_time);
            connes_apodization((Array){.data = dipy, .n = params->MaxTrajectoryLength}, params->sampling_time);
            connes_apodization((Array){.data = dipz, .n = params->MaxTrajectoryLength}, params->sampling_time);
            
            gsl_fft_real_radix2_transform(dipx, 1, params->MaxTrajectoryLength);
            gsl_fft_real_radix2_transform(dipy, 1, params->MaxTrajectoryLength);
            gsl_fft_real_radix2_transform(dipz, 1, params->MaxTrajectoryLength);

            gsl_fft_square(dipx, params->MaxTrajectoryLength);
            gsl_fft_square(dipy, params->MaxTrajectoryLength);
            gsl_fft_square(dipz, params->MaxTrajectoryLength);
               
            for (size_t i = 0; i < frequency_array_length; ++i) {
                sf_iter.data[i] += SF_COEFF * pr_mu * (dipx[i] + dipy[i] + dipz[i]);
            }

            memset(dipx, 0, params->MaxTrajectoryLength * sizeof(double)); 
            memset(dipy, 0, params->MaxTrajectoryLength * sizeof(double)); 
            memset(dipz, 0, params->MaxTrajectoryLength * sizeof(double)); 

            traj_counter++;
        } 

        MPI_Allreduce(MPI_IN_PLACE, sf_iter.data, frequency_array_length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        for (size_t i = 0; i < frequency_array_length; ++i) {
            sf_total.data[i] += sf_iter.data[i];
        }
        sf_total.ntraj += local_ntrajectories * _wsize;

        double M0_est = compute_M0_from_sf(sf_total) / sf_total.ntraj;
        double M2_est = compute_M2_from_sf(sf_total) / sf_total.ntraj;

        PRINT0("ITERATION %zu/%zu: accumulated %zu trajectories. Saving temporary result to '%s'\n", iter+1, params->niterations, sf_total.ntraj, params->sf_filename);
        PRINT0("M0 ESTIMATE FROM SF: %.5e, PRELIMINARY M0 ESTIMATE: %.5e, diff: %.3f%%\n",   M0_est, prelim_M0, (M0_est - prelim_M0)/prelim_M0*100.0);
        PRINT0("M2 ESTIMATE FROM SF: %.5e, PRELIMINARY M2 ESTIMATE: %.5e, diff: %.3f%%\n\n", M2_est, prelim_M2, (M2_est - prelim_M2)/prelim_M2*100.0);

        if (_wrank == 0) {
            save_spectral_function(fp, sf_total, params);
        }
    }

    free_sfnc(sf_iter);

    for (size_t i = 0; i < frequency_array_length; ++i) {
        sf_total.data[i] /= sf_total.ntraj;
    }

    return sf_total; 

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

void save_correlation_function(FILE *fp, CFnc cf, CalcParams *params)
/*
 * Here we assume that values of the correlation function come in unnormalized by the number of trajectories
 * (which is set in clrn.ntraj)
 */ 
{
    assert(cf.ntraj > 0);

    // truncate the file to a length of 0, effectively clearing its contents
    int fd = fileno(fp); 
    if (ftruncate(fd, 0) < 0) {
        printf("ERROR: could not truncate file: %s\n", strerror(errno));
        exit(1);
    }

    // resets the file position indicator
    rewind(fp);

    fprintf(fp, "# HAWAII HYBRID v0.1\n");
    fprintf(fp, "# TEMPERATURE: %.2f\n", cf.Temperature);
    fprintf(fp, "# PAIR STATE: %s\n", pair_state_name(params->ps));
    fprintf(fp, "# AVERAGE OVER %zu TRAJECTORIES\n", cf.ntraj); 
    fprintf(fp, "# MAXIMUM TRAJECTORY LENGTH: %zu\n", cf.len);
    fprintf(fp, "# PARTIAL PARTITION FUNCTION: %.8e\n", params->partial_partition_function_ratio);
    fprintf(fp, "# CVODE TOLERANCE: %.2e\n", params->cvode_tolerance);

    for (size_t i = 0; i < cf.len; ++i) {
        fprintf(fp, "%.2f %.10e\n", cf.t[i], cf.data[i] / cf.ntraj * ALU*ALU*ALU);
    }
    
    fflush(fp);
}
            
void save_spectral_function(FILE *fp, SFnc sf, CalcParams *params) 
/*
 * Here we assume that values of the spectral function come in unnormalized by the number of trajectories
 * (which is set in sf.ntraj)
 */ 
{
    assert(sf.ntraj > 0);

    // truncate the file to a length of 0, effectively clearing its contents
    int fd = fileno(fp); 
    if (ftruncate(fd, 0) < 0) {
        printf("ERROR: could not truncate file: %s\n", strerror(errno));
        exit(1);
    }
    
    // resets the file position indicator
    rewind(fp);

    fprintf(fp, "# HAWAII HYBRID v0.1\n");
    fprintf(fp, "# TEMPERATURE: %.2f\n", sf.Temperature);
    fprintf(fp, "# AVERAGE OVER %zu TRAJECTORIES\n", sf.ntraj); 
    fprintf(fp, "# MAXIMUM TRAJECTORY LENGTH: %zu\n", params->MaxTrajectoryLength);
    fprintf(fp, "# INITIAL DISTANCE: %.3f\n", params->R0);
    fprintf(fp, "# CVODE TOLERANCE: %.2e\n", params->cvode_tolerance);

    for (size_t i = 0; i < sf.len; ++i) {
        fprintf(fp, "%.10f   %.10e\n", sf.nu[i], sf.data[i] / sf.ntraj); 
    }

    fflush(fp);
}


/*
#define nob_da_append_many(da, new_items, new_items_count)                                  \
    do {                                                                                    \
        if ((da)->count + new_items_count > (da)->capacity) {                               \
            if ((da)->capacity == 0) {                                                      \
                (da)->capacity = NOB_DA_INIT_CAP;                                           \
            }                                                                               \
            while ((da)->count + new_items_count > (da)->capacity) {                        \
                (da)->capacity *= 2;                                                        \
            }                                                                               \
            (da)->items = NOB_REALLOC((da)->items, (da)->capacity*sizeof(*(da)->items));    \
            NOB_ASSERT((da)->items != NULL && "ASSERT: not enough memory!\n");              \
        }                                                                                   \
        memcpy((da)->items + (da)->count, new_items, new_items_count*sizeof(*(da)->items)); \
        (da)->count += new_items_count;                                                     \
    } while (0)
*/

void sb_append(String_Builder *sb, const char *line, size_t n) {
    strncpy(sb->items + sb->count, line, n);
    sb->count += n;
}

void free_sb(String_Builder sb) {
    free(sb.items);
}

bool read_correlation_function(const char *filename, String_Builder *sb, CFnc *cf) 
{
    bool result = true;

    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: could not open the file '%s': %s\n", filename, strerror(errno));
        return_defer(false); 
    }
  
    char* line = NULL;
    size_t n;

    size_t header_lines = 0;

    while (getline(&line, &n, fp) > 0) {
        n = strlen(line);

        if (line[0] != '#') break;

        size_t new_count = sb->count + n;
        if (new_count > sb->capacity) {
            sb->items = (char*) realloc(sb->items, new_count);
            assert((sb->items != NULL) && "ASSERT: not enough memory!\n");
            sb->capacity = new_count;
        }

        sb_append(sb, line, n);
    
        free(line);
        line = NULL; 

        header_lines++;
    }
   
    printf("# of lines in header: %zu\n", header_lines);

    // rewind the file pointer by one line because we read forward while reading header 
    long pos = ftell(fp);
    fseek(fp, pos - strlen(line), SEEK_SET);

    size_t INIT_CAPACITY = 4096;

    cf->t    = (double*) malloc(INIT_CAPACITY * sizeof(double));
    cf->data = (double*) malloc(INIT_CAPACITY * sizeof(double));
    cf->len = 0;
    cf->capacity = INIT_CAPACITY;

    size_t lineno = header_lines + 1;
    double vt, vc;

    while (fscanf(fp, "%lf %lf", &vt, &vc) == 2) {
        if (cf->len + 1 > cf->capacity) {
            size_t new_capacity = 2 * cf->capacity;

            cf->t = (double*) realloc(cf->t, new_capacity * sizeof(double));
            assert((cf->t != NULL) && "ASSERT: not enough memory!");
            cf->data = (double*) realloc(cf->data, new_capacity * sizeof(double));
            assert((cf->data != NULL) && "ASSERT: not enough memory!");
            cf->capacity = new_capacity;
        }

        if (isnan(vt)) {
            printf("ERROR: could not parse time value on line %zu!\n", lineno);
            return_defer(false); 
        }

        if (isnan(vc)) {
            printf("ERROR: could not parse correlation function value on line %zu!\n", lineno);
            return_defer(false); 
        }
     
        cf->t[cf->len]    = vt;
        cf->data[cf->len] = vc;
        cf->len++;
        lineno++;
    }

defer:
    if (fp) fclose(fp);
    return result;
}

bool read_spectral_function(const char *filename, String_Builder *sb, SFnc *sf) 
{
    bool result = true;

    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: could not open the file '%s': %s\n", filename, strerror(errno));
        return_defer(false); 
    }
  
    char* line = NULL;
    size_t n;

    size_t header_lines = 0;

    while (getline(&line, &n, fp) > 0) {
        n = strlen(line);

        if (line[0] != '#') break;

        size_t new_count = sb->count + n;
        if (new_count > sb->capacity) {
            sb->items = (char*) realloc(sb->items, new_count);
            assert((sb->items != NULL) && "ASSERT: not enough memory!\n");
            sb->capacity = new_count;
        }

        sb_append(sb, line, n);

        free(line);
        line = NULL; 

        header_lines++;
    }
   
    printf("# of lines in header: %zu\n", header_lines);

    // rewind the file pointer by one line because we read forward while reading header 
    long pos = ftell(fp);
    fseek(fp, pos - strlen(line), SEEK_SET);

    size_t INIT_CAPACITY = 4096;

    sf->nu   = (double*) malloc(INIT_CAPACITY * sizeof(double));
    sf->data = (double*) malloc(INIT_CAPACITY * sizeof(double));
    sf->len = 0;
    sf->capacity = INIT_CAPACITY;

    size_t lineno = header_lines + 1;
    double vnu, vsf;

    while (fscanf(fp, "%lf %lf", &vnu, &vsf) == 2) {
        if (sf->len + 1 > sf->capacity) {
            size_t new_capacity = 2 * sf->capacity;

            sf->nu = (double*) realloc(sf->nu, new_capacity * sizeof(double));
            assert((sf->nu != NULL) && "ASSERT: not enough memory!");
            sf->data = (double*) realloc(sf->data, new_capacity * sizeof(double));
            assert((sf->data != NULL) && "ASSERT: not enough memory!");
            sf->capacity = new_capacity;
        }

        if (isnan(vnu)) {
            printf("ERROR: could not parse frequency value on line %zu!\n", lineno);
            return_defer(false); 
        }

        if (isnan(vsf)) {
            printf("ERROR: could not parse spectral function value on line %zu!\n", lineno);
            return_defer(false); 
        }
     
        sf->nu[sf->len] = vnu;
        sf->data[sf->len] = vsf;
        sf->len++;
        lineno++;
    }

defer:
    if (fp) fclose(fp);
    return result;
}

bool writetxt(const char *filename, double *x, double *y, size_t len, const char *header) 
{
    bool result = true;

    FILE *fp = fopen(filename, "wb");
    if (fp == NULL) {
        printf("ERROR: could not open the file '%s': %s\n", filename, strerror(errno));
        return_defer(false); 
    }

    size_t nchars = 0;
    nchars += fwrite(header, strlen(header), 1, fp);

    for (size_t i = 0; i < len; ++i) {
        nchars += fprintf(fp, "%.3f %.10e\n", x[i], y[i]);
    }

    printf("INFO: writing %zu characters into '%s' finished.\n", nchars, filename); 

defer:
    if (fp) fclose(fp);
    return result;
}

double analytic_full_partition_function_by_V(MoleculeSystem *ms, double Temperature)
{
    double pf_analytic = 0.0;

    if ((ms->m1.t == ATOM) && (ms->m2.t == ATOM)) {
        TODO("analytic_full_partition_function_by_V");  
    } else if ((ms->m1.t == LINEAR_MOLECULE) && (ms->m2.t == ATOM)) {
        pf_analytic = 4.0*M_PI * pow(2.0*M_PI*Temperature/HkT, 2.5) * pow(ms->mu, 1.5) * ms->m1.II[0];
    } else if ((ms->m1.t == ROTOR) && (ms->m2.t == LINEAR_MOLECULE)) {
        pf_analytic = 32.0*M_PI*M_PI*M_PI * pow(2.0*M_PI*Temperature/HkT, 4.0) * pow(ms->mu, 1.5) * sqrt(ms->m1.II[0]*ms->m1.II[1]*ms->m1.II[2]) * ms->m2.II[0];
    } else {
        TODO("analytic_full_partition_function_by_V");  
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

void free_sfnc(SFnc sf) {
    free(sf.nu);
    free(sf.data);
}

void free_spectrum(Spectrum sp) {
    free(sp.nu);
    free(sp.data);
}

WingParams INIT_WP = {
  .A = 1.0, 
  .B = 2.0, 
  .C = 3.0,
};

double wingmodel(WingParams *wp, double t) {
    return wp->C + wp->A / (1.0 + wp->B * wp->B * t * t);
}

int wingmodel_f(const gsl_vector* x, void* data, gsl_vector* f)
/*
 * [input]  x    : parameters
 * [input]  data : CF
 * [output] f    : Yi - yi = difference between model and CF values
 *
 * MODEL: Lorentzian function shifted upwards by constant: f = c + a /(1 + b^2 x^2)
 * https://mathworld.wolfram.com/LorentzianFunction.html 
 */ 
{
    WingData *wd = (WingData*) data; 
    size_t n  = wd->n;
    double* t = wd->t;
    double* y = wd->y;

    double A = gsl_vector_get(x, 0);
    double B = gsl_vector_get(x, 1);
    double C = gsl_vector_get(x, 2);

    for (size_t k = 0; k < n; ++k) {
        double Yk = C + A / (1.0 + B * B * t[k] * t[k]);
        gsl_vector_set(f, k, Yk - y[k]);
    }

    return GSL_SUCCESS;
}

int wingmodel_df(const gsl_vector* x, void* data, gsl_matrix * J)
/*
 * [input]  x    : parameters
 * [input]  data : CF
 * [output] J    : Jacobian matrix J(i, j) = dfi / dxj
 *                 where fi = (Yi - yi) -- difference between model and CF values
 *                       and the xj are the parameters  
 *
 * MODEL: Lorentzian function shifted upwards by constant: f = c + a /(1 + b^2 x^2)
 * https://mathworld.wolfram.com/LorentzianFunction.html 
 */ 
{
    size_t n = ((WingData*) data)->n;
    double* t = ((WingData*) data)->t;

    double A = gsl_vector_get(x, 0);
    double B = gsl_vector_get(x, 1);

    for (size_t k = 0; k < n; ++k) {
        gsl_matrix_set(J, k, 0,  1.0 / (1.0 + B * B * t[k] * t[k]));
        gsl_matrix_set(J, k, 1, -2.0 * A * B * t[k] * t[k] / (B * B * t[k] * t[k] + 1) / (B * B * t[k] * t[k] + 1));
        gsl_matrix_set(J, k, 2,  1.0);
    }

    return GSL_SUCCESS;
}

void gsl_multifit_callback(const size_t iter, void* params, const gsl_multifit_nlinear_workspace* w)
{
    UNUSED(params);

    gsl_vector* f = gsl_multifit_nlinear_residual(w);
    gsl_vector* x = gsl_multifit_nlinear_position(w);
    double rcond;

    // compute reciprocal condition number of J(x) 
    gsl_multifit_nlinear_rcond(&rcond, w);

    fprintf(stdout, "    [callback] LM iter %2zu: A = %.12f, B = %.12f, C = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n",
            iter, gsl_vector_get(x, 0), gsl_vector_get(x, 1), gsl_vector_get(x, 2), 1.0 / rcond, gsl_blas_dnrm2(f));
}

void gsl_nonlinear_opt(size_t n, double* x, double* y, WingParams *wing_params)
/*
 * fit wing model to [x, y] data using the Levenberg-Marquardt method
 *
 * based on example from
 * https://www.gnu.org/software/gsl/doc/html/nls.html 
 *
 * WingParams = {A, B, C} : initial parameters provided to LM method
 *                          optimized parameters are written back into array 
 *
 */
{
#define nparams 3 

    printf("n = %zu\n", n);

    WingData d = {
        .n = n, 
        .t = x, 
        .y = y,
    };

    const gsl_multifit_nlinear_type* T = gsl_multifit_nlinear_trust;

    gsl_multifit_nlinear_fdf fdf;
    fdf.f      = wingmodel_f;
    fdf.df     = wingmodel_df; // set to NULL for finite-difference Jacobian
    fdf.fvv    = NULL; // not using geodesic acceleration
    fdf.n      = n;
    fdf.p      = nparams;
    fdf.params = (void*) &d;
    
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
    gsl_multifit_nlinear_workspace* w = gsl_multifit_nlinear_alloc(T, &fdf_params, n, nparams);
    
    double *weights = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; ++i) { 
        weights[i] = 1.0;
    }

    double _wp[nparams] = {wing_params->A, wing_params->B, wing_params->C};
    gsl_vector_view wp = gsl_vector_view_array(_wp, nparams);

    gsl_vector_view wts = gsl_vector_view_array(weights, n);
    gsl_multifit_nlinear_winit(&wp.vector, &wts.vector, &fdf, w); // initialize solver with starting point and weights

    // compute initial cost function
    gsl_vector* f = gsl_multifit_nlinear_residual(w);

    double chisq0;
    gsl_blas_ddot(f, f, &chisq0);

    int status, info;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 1e-8;

    const size_t niter_max = 100;
    status = gsl_multifit_nlinear_driver(niter_max, xtol, gtol, ftol, gsl_multifit_callback, NULL, &info, w);

    // covariance of best fit parameters
    gsl_matrix* J = gsl_multifit_nlinear_jac(w);
    gsl_matrix* covar = gsl_matrix_alloc(3, 3);
    gsl_multifit_nlinear_covar(J, 0.0, covar);

    double chisq;
    gsl_blas_ddot(f, f, &chisq); // final cost

    fprintf(stderr, "    summary from method '%s/%s'\n", gsl_multifit_nlinear_name(w), gsl_multifit_nlinear_trs_name(w));
    fprintf(stderr, "    number of iterations: %zu\n", gsl_multifit_nlinear_niter(w));
    fprintf(stderr, "    function evaluations: %zu\n", fdf.nevalf);
    fprintf(stderr, "    Jacobian evaluations: %zu\n", fdf.nevaldf);
    fprintf(stderr, "    reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "    initial |f(x)| = %f\n", sqrt(chisq0));
    fprintf(stderr, "    final   |f(x)| = %f\n", sqrt(chisq));

    double dof = n - nparams;
    double c = fmax(1.0, sqrt(chisq / dof));
    
    wing_params->A = gsl_vector_get(w->x, 0);
    wing_params->B = gsl_vector_get(w->x, 1);
    wing_params->C = gsl_vector_get(w->x, 2);

    fprintf(stderr, "chisq / dof = %g\n", chisq / dof);
    fprintf(stderr, "A = %.10f +/- %.10f\n", wing_params->A, c * sqrt(gsl_matrix_get(covar, 0, 0)));
    fprintf(stderr, "B = %.10f +/- %.10f\n", wing_params->B, c * sqrt(gsl_matrix_get(covar, 1, 1)));
    fprintf(stderr, "C = %.10f +/- %.10f\n", wing_params->C, c * sqrt(gsl_matrix_get(covar, 2, 2)));

    fprintf(stderr, "status = %s\n", gsl_strerror(status));

    gsl_multifit_nlinear_free(w);
    gsl_matrix_free(covar);

    free(weights);

#undef nparams 
}

WingParams fit_baseline(CFnc *cf, size_t EXT_RANGE_MIN)
// EXT_RANGE_MIN :: leftmost value of `tail` range to fit the parameters of the baseline function

// See paper "Simulation of collision-induced absorption spectra based on classical trajectories and
// ab initio potential and induced dipole surfaces. II. CO2-Ar rototranslational band including true 
// dimer contribution"
//
// we rely on the specific structure of the correlation function (see Fig. 7)
// specifically, we expect it to have a secondary maximum from which it decays to some POSITIVE constant
{
    double CFmin  = FLT_MAX; 
    double CFmax  = FLT_MIN; 

    // auxiliary variable to detect that we traversed the first minimum of the CF
    bool wasnegative = false;

    // secondary maximum: its time coordinate [max2time] and CF value [CFmax2]
    double max2time = 0.0; 
    double CFmax2 = FLT_MIN; 

    for (size_t k = 0; k < cf->len; ++k) {
        if (cf->data[k] > CFmax) CFmax = cf->data[k];
        if (cf->data[k] < CFmin) CFmin = cf->data[k];

        if (cf->data[k] < 0.0) wasnegative = true;
        if (wasnegative && (cf->data[k] > CFmax2)) { CFmax2 = cf->data[k]; max2time = cf->t[k]; }
    } 

    printf("INFO: CFmin    = %.10e\n", CFmin); 
    printf("INFO: CFmax    = %.10e\n", CFmax);
    printf("INFO: max2time = %.10e\n", max2time / ATU); 
    printf("INFO: CFmax2   = %.10e\n", CFmax2);
    
    assert(cf->len > EXT_RANGE_MIN);
    
    size_t itotal = cf->len - EXT_RANGE_MIN + 1;
    double *xdat2 = (double*) malloc(itotal * sizeof(double));
    double *ydat2 = (double*) malloc(itotal * sizeof(double));
   
    // selecting values of `wing tail` and also scaling them  
    for (size_t k = EXT_RANGE_MIN, ind = 0; k < cf->len; ++k, ++ind) {
        xdat2[ind] = cf->t[k] / max2time;
        ydat2[ind] = cf->data[k] / CFmax2;
    }

    WingParams wp = {
        .A = INIT_WP.A, 
        .B = INIT_WP.B, 
        .C = fmax(ydat2[itotal-1], INIT_WP.C)
    }; 
    
    printf("Initial parameter values: %.5e %.5e %.5e\n", wp.A, wp.B, wp.C);

    gsl_nonlinear_opt(itotal, xdat2, ydat2, &wp);

    //assert(wp.A > 0);
    //assert(wp.B > 0);
    //assert(wp.C > 0);
    if (wp.C < 0) {
        printf("\n");
        printf("!!! WingParams.C is negative! Deal with it. Continuing...\n\n");
    }

    if (wp.B < 0) {
        printf("\n");
        printf("!!! WingParams.B is negative! Probably the Lorentzian model is not applicable in selected range. Consider fitting only the constant. Continuing...\n\n");
    }

    max2time = max2time / ATU;

    // descale parameters after optimization
    wp.A *= CFmax2;
    wp.B /= max2time;
    wp.C *= CFmax2;

    free(xdat2);
    free(ydat2);

    return wp; 
}

void connes_apodization(Array a, double sampling_time) 
{
    double tmax = sampling_time * a.n;

    for (size_t i = 0; i < a.n; ++i) {
        double tcurr = sampling_time * i;
        a.data[i] *= (1.0 - (tcurr / tmax) * (tcurr / tmax)) * (1.0 - (tcurr / tmax) * (tcurr / tmax)); 
    } 
}

#define REAL(z, i) ((z)[2*(i)])
#define IMAG(z, i) ((z)[2*(i) + 1])

double *idct(double *v, size_t len)
/*
 * Implementation of Inverse Fast Discrete Cosine Transform using the same Makhoul trick (2N + half-sample shift) and IFFT.
 * See original paper:
 * http://eelinux.ee.usm.maine.edu/courses/ele486/docs/makhoul.fastDCT.pdf 
 * 
 * The packing of input array for IFFT procedure is more clearly demonstrated by the following Python code: 
 * import numpy as np
 * import scipy as sp

 * N = 4
 * v = np.array([1.0, 2.0, 3.0, 4.0])
 * f = sp.fft.dct(v)
 *
 * f_mkh = np.zeros(2*N).astype(complex)
 * f_mkh[:N] = f 
 * f_mkh[N] = 0.0
 * f_mkh[N + 1:] = -f[1:][::-1]

 * ce = np.array([np.exp(1j * np.pi * k / (2.0 * N)) for k in range(2 * N)])
 * f_mkh = f_mkh * ce
 * v_mkh = np.fft.ifft(f_mkh).real[:N]
 * print(np.allclose(v, v_mkh))
 */
{
    assert(is_power_of_two(len) && "Length of the input array should be a power of 2"); 
  
    double *packed_v_mkh = (double*) malloc(4*len * sizeof(double)); 
    memset(packed_v_mkh, 0, 4*len*sizeof(double));

    double complex cx;

    for (size_t k = 0; k < len; ++k) {
        cx = v[k] * cexp(I * M_PI * (double)k / (2.0 * len)) / (2.0 * len);

        if (isnan(creal(cx)) || isnan(cimag(cx))) {
            printf("ERROR: idct: the value calculated at k = %zu is NaN! Check the provided array\n", k);
            exit(1);
        }

        REAL(packed_v_mkh, k) = creal(cx); 
        IMAG(packed_v_mkh, k) = cimag(cx);
    }
     
    for (size_t k = len + 1; k < 2 * len; ++k) {
        size_t i = 2*len - k;
        cx = -v[i] * cexp(I * M_PI * (double)k / (2.0 * len)) / (2.0 * len); 
        REAL(packed_v_mkh, k) = creal(cx); 
        IMAG(packed_v_mkh, k) = cimag(cx); 
        
        if (isnan(creal(cx)) || isnan(cimag(cx))) {
            printf("ERROR: idct: the value calculated at k = %zu is NaN! Check the provided array\n", k);
            exit(1);
        }
    }

    // Note:
    // the length of array passed to gsl_fft_complex_radix2... is equal to number of `complex numbers` in it = 2 * sz
    // even though there are 4 * sz `double numbers` 
    gsl_fft_complex_radix2_backward(packed_v_mkh, 1, 2*len);
    
    double *f = (double*) malloc(len * sizeof(double));
    memset(f, 0, len*sizeof(double));

    for (size_t k = 0; k < len; ++k) {
        f[k] = REAL(packed_v_mkh, k);
    } 

    free(packed_v_mkh);

    return f;
}

SFnc dct_numeric_sf(CFnc cf, WingParams *wp)
{
    double Xscale = 1.0 / LightSpeed_cm / ATU / 2.0 / M_PI;
    double Yscale = ATU * ADIPMOMU * ADIPMOMU / (4.0 * M_PI * EPSILON0);

    double *cfnum = malloc(cf.len * sizeof(double));
    for (size_t i = 0; i < cf.len; ++i) {
        cfnum[i] = cf.data[i] - wingmodel(wp, cf.t[i]);    
    }
    
    double dt = (cf.t[1] - cf.t[0]) / ATU;
    connes_apodization((Array) {.data = cfnum, .n = cf.len }, dt);

    double *Ft = idct(cfnum, cf.len);
    free(cfnum);

    double tmax = cf.len * dt;

    SFnc sfnum = {
        .nu   = malloc(cf.len * sizeof(double)),
        .data = malloc(cf.len * sizeof(double)),
        .len  = cf.len,
    };

    for (size_t i = 0; i < cf.len; ++i) {
        sfnum.nu[i]   = i * M_PI / tmax * Xscale;
        sfnum.data[i] = Ft[i] * tmax / M_PI * Yscale;
    }
  
    free(Ft);

    return sfnum; 
}

SFnc desymmetrize_sch(SFnc sf) 
{
    SFnc sfd = {
        .nu   = malloc(sf.len * sizeof(double)),
        .data = malloc(sf.len * sizeof(double)),
        .len  = sf.len,
    };

    memcpy(sfd.nu, sf.nu, sf.len * sizeof(double));

    for (size_t i = 0; i < sf.len; ++i) {
        double hnukt = Planck * LightSpeed_cm * sf.nu[i] / (Boltzmann * sf.Temperature);
        sfd.data[i] = sf.data[i] * exp(hnukt / 2.0);
    }

    return sfd;
}

Spectrum compute_alpha(SFnc sf) 
{
    Spectrum sp = {
        .nu   = malloc(sf.len * sizeof(double)),
        .data = malloc(sf.len * sizeof(double)),
        .len  = sf.len,
    };

    memcpy(sp.nu, sf.nu, sf.len * sizeof(double));

    for (size_t i = 0; i < sf.len; ++i) {
        double hnukt = Planck * LightSpeed_cm * sp.nu[i] / (Boltzmann * sf.Temperature); 
        sp.data[i] = Moment_SF_Coeff * sp.nu[i] * (1.0 - exp(-hnukt)) * sf.data[i]; 
    }

    return sp;    
}

double integrate_composite_simpson(double *x, double *y, size_t len) 
/*
 * Composite Simpson's 3/8 rule
 */
{
    double h = x[1] - x[0];
    double sum = y[0] + y[len - 1]; 
    
    for (size_t j = 1; j < len - 1; ++j) {
        if (j % 3 == 0) {
            sum += 2.0 * y[j];
        } else {
            sum += 3.0 * y[j];
        }
    } 
    
    return sum * 3.0 * h / 8.0;
} 

double compute_M0_from_sf(SFnc sf) {
    return 2.0*Moment_SF_Coeff * integrate_composite_simpson(sf.nu, sf.data, sf.len);
}

double compute_M2_from_sf(SFnc sf) {
    double *y = (double*) malloc(sf.len * sizeof(double));

    for (size_t i = 0; i < sf.len; ++i) {
        y[i] = sf.nu[i] * sf.nu[i] * sf.data[i];
    }

    double M2 = 2.0*Moment_SF_Coeff * integrate_composite_simpson(sf.nu, y, sf.len);
    free(y);

    return M2;
}
