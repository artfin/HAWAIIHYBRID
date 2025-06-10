#include "hep_hawaii.hpp"
#include "hep_hawaii.h"

static MoleculeSystem* gms = NULL;
static CalcParams *gparams = NULL;
static double gT = 0.0;

void transform_variables(hep::mc_point<double> const& x, double* qp, double* jac)
/*
 * Input: n-dimensional point in [0; 1]^n  hypercube
 * Output: [std::vector<double>] point in phase space
 *         [double]              cumulative jacobian of the transform
 */
{
    assert(x.point().size() == gms->QP_SIZE);

    *jac = 1;

    qp[IR] = 1.0 / x.point()[0]; 
    *jac *= qp[IR] * qp[IR];  
    
    qp[IPR] =  std::tan(M_PI * (x.point()[IPR] - 0.5));
    *jac *= M_PI * (1.0 + qp[IPR] * qp[IPR]);
    
    qp[IPHI] = x.point()[IPHI] * 2.0 * M_PI;
    *jac *= 2.0 * M_PI;
    
    qp[IPPHI] =  std::tan(M_PI * (x.point()[IPPHI] - 0.5));
    *jac *= M_PI * (1.0 + qp[IPPHI] * qp[IPPHI]);
    
    qp[ITHETA] = x.point()[ITHETA] * M_PI;
    *jac *= M_PI;
    
    qp[IPTHETA] =  std::tan(M_PI * (x.point()[IPTHETA] - 0.5));
    *jac *= M_PI * (1.0 + qp[IPTHETA] * qp[IPTHETA]);

    switch (gms->m1.t) {
        case ATOM: break;
        case LINEAR_MOLECULE_REQ_HALFINTEGER:
        case LINEAR_MOLECULE_REQ_INTEGER:
        case LINEAR_MOLECULE: {
          qp[6 + IPHI] = x.point()[6 + IPHI] * 2.0 * M_PI;
          *jac *= 2.0 * M_PI;
        
          qp[6 + IPPHI] = std::tan(M_PI * (x.point()[6 + IPPHI] - 0.5));
          *jac *= M_PI * (1.0 + qp[6 + IPPHI] * qp[6 + IPPHI]);
          
          qp[6 + ITHETA] = x.point()[6 + ITHETA] * M_PI;
          *jac *= M_PI;
        
          qp[6 + IPTHETA] = std::tan(M_PI * (x.point()[6 + IPTHETA] - 0.5));
          *jac *= M_PI * (1.0 + qp[6 + IPTHETA] * qp[6 + IPTHETA]);
          
          break;
        }
        case ROTOR: {
          qp[6 + IPHI] = x.point()[6 + IPHI] * 2.0 * M_PI;
          *jac *= 2.0 * M_PI;
        
          qp[6 + IPPHI] = std::tan(M_PI * (x.point()[6 + IPPHI] - 0.5));
          *jac *= M_PI * (1.0 + qp[6 + IPPHI] * qp[6 + IPPHI]);
          
          qp[6 + ITHETA] = x.point()[6 + ITHETA] * M_PI;
          *jac *= M_PI;
        
          qp[6 + IPTHETA] = std::tan(M_PI * (x.point()[6 + IPTHETA] - 0.5));
          *jac *= M_PI * (1.0 + qp[6 + IPTHETA] * qp[6 + IPTHETA]);

          qp[6 + IPSI] = x.point()[6 + IPSI] * 2.0 * M_PI;
          *jac *= 2.0 * M_PI;
        
          qp[6 + IPPSI] = std::tan(M_PI * (x.point()[6 + IPPSI] - 0.5));
          *jac *= M_PI * (1.0 + qp[6 + IPPSI] * qp[6 + IPPSI]);
          
          break;
        }
        case ROTOR_REQUANTIZED_ROTATION: {
          assert(0 && "ERROR: not applicable for ROTOR_REQUANTIZED_ROTATION\n"); // ? 
        }
    }

    switch (gms->m2.t) {
        case ATOM: break; 
        case LINEAR_MOLECULE: {
          qp[6 + gms->m1.t%MODULO_BASE + IPHI] = x.point()[6 + gms->m1.t%MODULO_BASE + IPHI] * 2.0 * M_PI;
          *jac *= 2.0 * M_PI;
        
          qp[6 + gms->m1.t%MODULO_BASE + IPPHI] = std::tan(M_PI * (x.point()[6 + gms->m1.t%MODULO_BASE + IPPHI] - 0.5));
          *jac *= M_PI * (1.0 + qp[6 + gms->m1.t%MODULO_BASE + IPPHI] * qp[6 + gms->m1.t%MODULO_BASE + IPPHI]);
          
          qp[6 + gms->m1.t%MODULO_BASE + ITHETA] = x.point()[6 + gms->m1.t%MODULO_BASE + ITHETA] * M_PI;
          *jac *= M_PI;
        
          qp[6 + gms->m1.t%MODULO_BASE + IPTHETA] = std::tan(M_PI * (x.point()[6 + gms->m1.t%MODULO_BASE + IPTHETA] - 0.5));
          *jac *= M_PI * (1.0 + qp[6 + gms->m1.t%MODULO_BASE + IPTHETA] * qp[6 + gms->m1.t%MODULO_BASE + IPTHETA]);

          break;
        }
        case ROTOR: {
          TODO("transform_variables");
        }
        case LINEAR_MOLECULE_REQ_INTEGER: {
          assert(0 && "ERROR: not applicable for LINEAR_MOLECULE_REQ_INTEGER\n"); // ? 
        }
        case LINEAR_MOLECULE_REQ_HALFINTEGER: {
          assert(0 && "ERROR: not applicable for LINEAR_MOLECULE_REQ_HALFINTEGER\n"); // ? 
        }
        case ROTOR_REQUANTIZED_ROTATION: {
          assert(0 && "ERROR: not applicable for ROTOR_REQUANTIZED_ROTATION\n"); // ? 
        }
    }
}

double integrand_pf(hep::mc_point<double> const& x)
{
    assert(gms != NULL);
    assert(gparams != NULL);
    assert(gT > 0);
    
    double jac;
    double qp[gms->QP_SIZE];
    transform_variables(x, qp, &jac);
    put_qp_into_ms(gms, (Array) {.data = qp, .n = gms->QP_SIZE});
    
    double R = qp[IR];
    if ((R < gparams->sampler_Rmin) || (R > gparams->sampler_Rmax)) {
        return 0.0;
    }
    
    double energy = Hamiltonian(gms); 
    
    if (gparams->ps == PAIR_STATE_FREE_AND_METASTABLE) {
        if (energy < 0.0) return 0.0; 
    }

    if (gparams->ps == PAIR_STATE_BOUND) {
        if (energy > 0.0) return 0.0; 
    }
    
    return jac * std::exp(-energy * HkT / gT);
}

double integrand_M0(hep::mc_point<double> const& x) 
{
    assert(dipole != NULL);
    assert(gms != NULL);
    assert(gparams != NULL);
    assert(gT > 0);

    double jac;
    double qp[gms->QP_SIZE];
    transform_variables(x, qp, &jac);
    put_qp_into_ms(gms, (Array) {.data = qp, .n = gms->QP_SIZE});

    double R = qp[IR];
    if (R < gparams->sampler_Rmin || R > gparams->sampler_Rmax) {
        return 0.0;
    }

    double energy = Hamiltonian(gms); 
    
    if (gparams->ps == PAIR_STATE_FREE_AND_METASTABLE) {
        if (energy < 0.0) return 0.0; 
    }

    if (gparams->ps == PAIR_STATE_BOUND) {
        if (energy > 0.0) return 0.0; 
    }
    
    double dip[3];
    extract_q_and_write_into_ms(gms);
    (*dipole)(gms->intermediate_q, dip);
    
    double dipsq = dip[0]*dip[0] + dip[1]*dip[1] + dip[2]*dip[2]; 

    //std::cout << "R: " << R << " => jac = " << jac << ", energy = " << energy << "\n";

    return jac * dipsq * std::exp(-energy * HkT / gT);
}

double integrand_M2(hep::mc_point<double> const& x)
{
    assert(dipole != NULL);
    assert(gms != NULL);
    assert(gparams != NULL);
    assert(gT > 0);
    
    double h = 1.0e-3;
    
    // TODO: move the memory allocation outside this function
    gsl_matrix *D           = gsl_matrix_alloc(gms->Q_SIZE, 3);
    gsl_matrix *dHdp        = gsl_matrix_alloc(1, gms->Q_SIZE);
    gsl_matrix *dip_lab_dot = gsl_matrix_alloc(1, 3);
    
    double jac;
    double qp[gms->QP_SIZE];
    transform_variables(x, qp, &jac);
    put_qp_into_ms(gms, (Array) {.data = qp, .n = gms->QP_SIZE});

    double R = qp[IR];
    if (R < gparams->sampler_Rmin || R > gparams->sampler_Rmax) {
        return 0.0;
    }

    double energy = Hamiltonian(gms); 
    
    if (gparams->ps == PAIR_STATE_FREE_AND_METASTABLE) {
        if (energy < 0.0) return 0.0; 
    }

    if (gparams->ps == PAIR_STATE_BOUND) {
        if (energy > 0.0) return 0.0; 
    }
            
    extract_q_and_write_into_ms(gms);
            
    double dp[3];
    double dm[3];

    for (size_t i = 0; i < gms->Q_SIZE; ++i) {
        double tmp = gms->intermediate_q[i];
        
        gms->intermediate_q[i] = tmp + h;
        (*dipole)(gms->intermediate_q, dp);
        
        gms->intermediate_q[i] = tmp - h;
        (*dipole)(gms->intermediate_q, dm);

        gsl_matrix_set(D, i, 0, (dp[0] - dm[0])/(2.0*h)); 
        gsl_matrix_set(D, i, 1, (dp[1] - dm[1])/(2.0*h)); 
        gsl_matrix_set(D, i, 2, (dp[2] - dm[2])/(2.0*h)); 

        gms->intermediate_q[i] = tmp;
    } 

    compute_dHdp(gms, dHdp); 

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, dHdp, D, 0.0, dip_lab_dot);
            
    double dipsq = gsl_matrix_get(dip_lab_dot, 0, 0)*gsl_matrix_get(dip_lab_dot, 0, 0) + \
                   gsl_matrix_get(dip_lab_dot, 0, 1)*gsl_matrix_get(dip_lab_dot, 0, 1) + \
                   gsl_matrix_get(dip_lab_dot, 0, 2)*gsl_matrix_get(dip_lab_dot, 0, 2);

    //std::cout << "R: " << R << " => jac = " << jac << ", energy = " << energy << "\n";
    
    gsl_matrix_free(D);
    gsl_matrix_free(dHdp);
    gsl_matrix_free(dip_lab_dot);

    return jac * dipsq * std::exp(-energy * HkT / gT);
}

void mpi_perform_integration(MoleculeSystem *ms, Integrand integrand, CalcParams *params, double Temperature, size_t niterations, size_t npoints, double *m, double *q)
{
    gms = ms;
    gparams = params;
    gT = Temperature;

    auto results = hep::mpi_vegas(
        MPI_COMM_WORLD,
        hep::make_integrand<double>(integrand, ms->QP_SIZE),
        std::vector<size_t>(niterations, npoints)
    );
    
    hep::mc_result<double> result = hep::accumulate<hep::weighted_with_variance>(results.begin() + 1, results.end());    
    // double chi_square_dof = hep::chi_square_dof<hep::weighted_with_variance>(results.begin() + 1, results.end());

    *m = result.value();
    *q = result.error();
}

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

void c_mpi_perform_integration(MoleculeSystem *ms, IntegrandType integrand_type, CalcParams *params, double Temperature, size_t niterations, size_t npoints, double *m, double *q)
{
    Integrand integrand = {};
    
    switch (integrand_type) {
        case INTEGRAND_M0: integrand = integrand_M0; break;
        case INTEGRAND_M2: integrand = integrand_M2; break;
        case INTEGRAND_PF: integrand = integrand_pf; break;
        default: UNREACHABLE("Unknown integrand type"); 
    }
    
    hep::mpi_vegas_callback<double>(hep::mpi_vegas_verbose_callback<double>);
    mpi_perform_integration(ms, integrand, params, Temperature, niterations, npoints, m, q);
}

#ifdef __cplusplus
}
#endif // __cplusplus

